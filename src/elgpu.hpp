#pragma once

// EL point particles on the device: nekRS's lpm_t (GPU findpts/interpolation,
// integration, migration) plus the two pieces it lacks, written as OKL
// kernels (src/elgpu.okl):
//
//   projection  transpose interpolation onto the GLL nodes of the particle's
//               element, gather-scatter, division by the assembled mass matrix
//               (an L2 projection, conservative to round-off), then a Gaussian
//               filter applied as explicit diffusion to pseudo-time sigma^2/2
//               with nekRS's weak Laplacian (conservative, needs no ghost
//               particles because the projection is already assembled);
//   coupling    FU = F / (1 - phi_p) to fluid->o_explicitTerms(), optionally
//               div u = -(1/phi_f) D phi_f/Dt to fluid->o_div.
//
// Particles are migrated every step so that each lives on the rank owning its
// element; interpolation and projection are then local. Compared with the
// ppiclF path (src/elibm.hpp) the only host traffic is the migration and the
// element/reference-coordinate list from findpts (two small arrays).
//
// Not here yet: particle-particle contact (needs a device cell list), IBM
// bodies, walls (particles leaving a non-periodic domain are deleted).

#include <algorithm>
#include <array>
#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "lpm.hpp"
#include "opSEM.hpp"

namespace elgpu
{

using vec3 = std::array<double, 3>;

struct Config {
  double mu = 0;
  double rhof = 0;
  vec3 grav{};
  std::array<bool, 3> periodic{};
  double filterWidth = 0; // Gaussian FWHM
  bool twoWay = true; // false: particles feel the fluid, the fluid not the particles
  bool qtl = true;
  bool qtlFromFlux = true; // d phi/dt from the projected volume flux, else finite difference
  bool zeroMeanDiv = true;
  // lpm_t's RK3 and RK4 evaluate the intermediate stages with the start-of-
  // step state (nekRS v26 lpm.cpp, integrateRK3/4: o_y_ instead of
  // o_ytmp_), which makes them first order for a velocity-dependent RHS;
  // RK2 is implemented correctly. Verified in cases/spray_jet2d (single
  // droplet vs the drag ODE).
  int rkOrder = 2;
};

struct Totals {
  vec3 fluidMomentum{};
  vec3 elMomentum{};
  vec3 elForceOnFluid{};
  vec3 elForceProjected{}; // after the filter
  vec3 source{};           // int FU dV
  double divRms = 0, divErrRms = 0;
  double elMass = 0;
  long long nEl = 0;
};

// Call from UDF_LoadKernels.
inline void registerKernels(occa::properties &kernelInfo)
{
  std::string f = __FILE__;
  platform->kernelRequests.add("elgpu", f.substr(0, f.rfind('/')) + "/elgpu.okl", kernelInfo);
}

class Coupling
{
public:
  Coupling(nrs_t *nrs, const Config &cfg) : nrs(nrs), cfg(cfg)
  {
    auto mesh = nrs->meshV;
    const auto offset = nrs->fluid->fieldOffset;

    rhsKernel = platform->kernelRequests.load("elgpu", "elgpuRHS");
    forceKernel = platform->kernelRequests.load("elgpu", "elgpuForce");
    wrapKernel = platform->kernelRequests.load("elgpu", "elgpuWrap");
    depositKernel = platform->kernelRequests.load("elgpu", "elgpuDeposit");
    sourceKernel = platform->kernelRequests.load("elgpu", "elgpuSource");
    weakLaplacianKernel = platform->kernelRequests.load("core-weakLaplacianHex3D");
    o_fieldOffsetScan = platform->device.malloc<dlong>(1);
    {
      const dlong zero = 0;
      o_fieldOffsetScan.copyFrom(&zero, 1);
    }

    std::vector<dfloat> x(mesh->Nlocal), y(mesh->Nlocal), z(mesh->Nlocal);
    mesh->o_x.copyTo(x.data(), mesh->Nlocal);
    mesh->o_y.copyTo(y.data(), mesh->Nlocal);
    mesh->o_z.copyTo(z.data(), mesh->Nlocal);
    const std::vector<dfloat> *c[3] = {&x, &y, &z};
    for (int k = 0; k < 3; k++) {
      lo[k] = *std::min_element(c[k]->begin(), c[k]->end());
      hi[k] = *std::max_element(c[k]->begin(), c[k]->end());
    }
    MPI_Allreduce(MPI_IN_PLACE, lo.data(), 3, MPI_DOUBLE, MPI_MIN, comm());
    MPI_Allreduce(MPI_IN_PLACE, hi.data(), 3, MPI_DOUBLE, MPI_MAX, comm());

    const double sigma = cfg.filterWidth / (2 * std::sqrt(2 * std::log(2.0)));
    const double tau = 0.5 * sigma * sigma; // heat kernel variance 2 tau = sigma^2

    o_one = platform->device.malloc<dfloat>(offset);
    platform->linAlg->fill(offset, 1.0, o_one);
    o_acc = platform->device.malloc<dfloat>(nAcc * offset);
    o_FU = platform->device.malloc<dfloat>(nrs->fluid->fieldOffsetSum);
    o_div = platform->device.malloc<dfloat>(offset);
    o_phifLast = platform->device.malloc<dfloat>(offset);
    platform->linAlg->fill(nrs->fluid->fieldOffsetSum, 0.0, o_FU);
    platform->linAlg->fill(offset, 0.0, o_div);
    platform->linAlg->fill(offset, 1.0, o_phifLast);
    volume = platform->linAlg->sum(mesh->Nlocal, mesh->o_LMM, comm());
    setupSmoothing(tau);

    nrs->userSource = [this](double) {
      this->nrs->fluid->o_explicitTerms().copyFrom(o_FU, this->nrs->fluid->fieldOffsetSum);
    };
    if (cfg.qtl) {
      nrs->userDivergence = [this](double) { this->nrs->fluid->o_div.copyFrom(o_div, this->nrs->meshV->Nlocal); };
    }

    particles = std::make_unique<lpm_t>(mesh, nrs->dt, 0.01, 5e-13);
    particles->setSolver("RK");
    particles->rkOrder(cfg.rkOrder);
    particles->setVerbosity(lpm_t::VerbosityLevel::None);
    particles->addInterpField("fluid_velocity", mesh->dim, nrs->fluid->fieldOffset, nrs->fluid->o_U);
    particles->addVariable(mesh->dim, "velocity");
    particles->addProp("dp");
    particles->addProp("rhop");
    particles->addProp("gid");
    particles->addProp("np"); // droplets per parcel
    particles->addUserData(this);
    particles->setUserRHS([](lpm_t *lpm, double, const deviceMemory<dfloat> &o_y, void *ud, deviceMemory<dfloat> &o_ydot) {
      static_cast<Coupling *>(ud)->rhs(lpm, o_y, o_ydot);
    });
  }

  const vec3 &domainLo() const { return lo; }
  const vec3 &domainHi() const { return hi; }
  int smoothingSteps() const { return nSmooth; }
  lpm_t &lpm() { return *particles; }

  // x, v: (3, n) column-major as in elibm; created on the calling rank.
  // weight = droplets per parcel (empty: 1).
  void initialize(const std::vector<double> &x,
                  const std::vector<double> &v,
                  const std::vector<double> &dp,
                  const std::vector<double> &rhop,
                  double t0,
                  const std::vector<double> &weight = {})
  {
    const int n = dp.size();
    std::vector<dfloat> y0, prop;
    pack(x, v, dp, rhop, weight, y0, prop);
    particles->initialize(n, t0, y0);
    if (n) {
      const int np = particles->nProps();
      for (int k = 0; k < np; k++) {
        auto name = std::vector<std::string>{"dp", "rhop", "gid", "np"}[k];
        particles->setProp(name, deviceMemory<dfloat>(std::vector<dfloat>(prop.begin() + k * n, prop.begin() + (k + 1) * n)));
      }
    }
    particles->migrate();
    updateSources(nullptr);
  }

  // Add parcels during the run; collective (every rank calls it, any count).
  void inject(const std::vector<double> &x,
              const std::vector<double> &v,
              const std::vector<double> &dp,
              const std::vector<double> &rhop,
              const std::vector<double> &weight = {})
  {
    std::vector<dfloat> y0, prop;
    pack(x, v, dp, rhop, weight, y0, prop);
    particles->addParticles(dp.size(), y0, prop);
  }

  void step(double time, int tstep)
  {
    if (tstep == 0) {
      return;
    }
    particles->integrate(time);
    wrap();
    particles->migrate();
    particles->deleteParticles();
    const double dt = nrs->dt[0];
    updateSources(&dt);
  }

  Totals totals()
  {
    auto mesh = nrs->meshV;
    auto &fluid = nrs->fluid;
    const auto offset = fluid->fieldOffset;
    Totals tot;
    for (int k = 0; k < 3; k++) {
      tot.fluidMomentum[k] =
          cfg.rhof * platform->linAlg->innerProd(mesh->Nlocal, mesh->o_LMM, fluid->o_U.slice(k * offset, mesh->Nlocal), comm());
      tot.source[k] = platform->linAlg->innerProd(mesh->Nlocal, mesh->o_LMM, o_FU.slice(k * offset, mesh->Nlocal), comm());
      tot.elForceProjected[k] = platform->linAlg->innerProd(mesh->Nlocal, mesh->o_LMM, o_acc.slice(k * offset, mesh->Nlocal), comm());
    }

    const int n = particles->numParticles();
    const int poff = particles->fieldOffset();
    std::vector<dfloat> v(3 * std::max(poff, 1)), q(nAcc * std::max(poff, 1)), rhop(std::max(poff, 1));
    if (n) {
      particles->getDOF("velocity").copyTo(v.data(), 3 * poff);
      o_q.copyTo(q.data(), nAcc * poff);
      particles->getProp("rhop").copyTo(rhop.data(), n);
    }
    std::vector<double> buf(8, 0.0);
    for (int i = 0; i < n; i++) {
      const double m = rhop[i] * q[i + 3 * poff];
      for (int k = 0; k < 3; k++) {
        buf[k] += m * v[i + k * poff];
        buf[3 + k] += q[i + k * poff];
      }
      buf[6] += m;
      buf[7] += 1;
    }
    MPI_Allreduce(MPI_IN_PLACE, buf.data(), buf.size(), MPI_DOUBLE, MPI_SUM, comm());
    for (int k = 0; k < 3; k++) {
      tot.elMomentum[k] = buf[k];
      tot.elForceOnFluid[k] = buf[3 + k];
    }
    tot.elMass = buf[6];
    tot.nEl = std::llround(buf[7]);

    if (cfg.qtl) {
      auto o_divU = opSEM::strongDivergence(mesh, offset, fluid->o_U);
      auto o_err = platform->deviceMemoryPool.reserve<dfloat>(mesh->Nlocal);
      platform->linAlg->axpbyz(mesh->Nlocal, 1.0, o_divU, -1.0, o_div, o_err);
      auto rms = [&](const occa::memory &o_a) {
        auto o_a2 = platform->deviceMemoryPool.reserve<dfloat>(mesh->Nlocal);
        platform->linAlg->axmyz(mesh->Nlocal, 1.0, o_a, o_a, o_a2);
        return std::sqrt(platform->linAlg->innerProd(mesh->Nlocal, mesh->o_LMM, o_a2, comm()) / volume);
      };
      tot.divRms = rms(o_div);
      tot.divErrRms = rms(o_err);
    }
    return tot;
  }

private:
  nrs_t *nrs;
  Config cfg;
  vec3 lo{}, hi{};
  double volume = 0;
  int nSmooth = 0;
  int lapSign = 0;
  double dtSmooth = 0;
  std::unique_ptr<lpm_t> particles;
  occa::kernel rhsKernel, forceKernel, wrapKernel, depositKernel, sourceKernel, weakLaplacianKernel;
  occa::memory o_fieldOffsetScan;
  static constexpr int nAcc = 7; // F(3), vol, vol*v(3)
  occa::memory o_one, o_acc, o_FU, o_div, o_phifLast;
  occa::memory o_q, o_el, o_r, o_Up; // per particle: (F, vol), element, reference coords, fluid velocity

  static MPI_Comm comm() { return platform->comm.mpiComm(); }

  long long nElGlobal = 0;

  // lpm_t layout: y0 (nDOF, n) and prop (nProps, n), n fastest
  void pack(const std::vector<double> &x,
            const std::vector<double> &v,
            const std::vector<double> &dp,
            const std::vector<double> &rhop,
            const std::vector<double> &weight,
            std::vector<dfloat> &y0,
            std::vector<dfloat> &prop)
  {
    const int n = dp.size();
    const int vId = particles->dofId("velocity");
    y0.assign(static_cast<size_t>(particles->nDOFs()) * n, 0.0);
    prop.assign(static_cast<size_t>(particles->nProps()) * n, 0.0);
    long long gid0 = n;
    MPI_Exscan(MPI_IN_PLACE, &gid0, 1, MPI_LONG_LONG, MPI_SUM, comm());
    if (platform->comm.mpiRank() == 0) {
      gid0 = 0;
    }
    for (int i = 0; i < n; i++) {
      for (int k = 0; k < 3; k++) {
        y0[i + k * n] = x[3 * i + k];
        y0[i + (vId + k) * n] = v[3 * i + k];
      }
      prop[i + particles->propId("dp") * n] = dp[i];
      prop[i + particles->propId("rhop") * n] = rhop[i];
      prop[i + particles->propId("gid") * n] = nElGlobal + gid0 + i;
      prop[i + particles->propId("np") * n] = weight.empty() ? 1.0 : weight[i];
    }
    long long total = n;
    MPI_Allreduce(MPI_IN_PLACE, &total, 1, MPI_LONG_LONG, MPI_SUM, comm());
    nElGlobal += total;
  }

  void rhs(lpm_t *lpm, const deviceMemory<dfloat> &o_y, deviceMemory<dfloat> &o_ydot)
  {
    lpm->interpolate();
    const int n = lpm->numParticles();
    if (n == 0) {
      return;
    }
    rhsKernel(n,
              lpm->fieldOffset(),
              lpm->dofId("velocity"),
              cfg.mu,
              cfg.rhof,
              cfg.grav[0],
              cfg.grav[1],
              cfg.grav[2],
              o_y,
              lpm->getProp("dp"),
              lpm->getProp("rhop"),
              lpm->getInterpField("fluid_velocity"),
              o_ydot);
  }

  void wrap()
  {
    const int n = particles->numParticles();
    if (n == 0) {
      return;
    }
    double len[3];
    for (int k = 0; k < 3; k++) {
      len[k] = cfg.periodic[k] ? hi[k] - lo[k] : -1.0;
    }
    wrapKernel(n, particles->fieldOffset(), lo[0], lo[1], lo[2], len[0], len[1], len[2], particles->getDOF("x"));
  }

  // L f = invM dssum(K f) with the SEM stiffness matrix K of the elliptic
  // solvers (the kernel behind opSEM::laplacian, which nekRS does not use and
  // which passes an uninitialised field-offset table to it). The weak form
  // integrates to zero over a periodic domain, so an explicit step
  // f += dt L f conserves int f dV.
  void applyLaplacian(const occa::memory &o_f, occa::memory &o_Lf)
  {
    auto mesh = nrs->meshV;
    weakLaplacianKernel(mesh->Nelements, 1, o_fieldOffsetScan, mesh->o_ggeo, mesh->o_D, o_one, o_f, o_Lf);
    oogs::startFinish(o_Lf, 1, 0, ogsDfloat, ogsAdd, mesh->oogs);
    platform->linAlg->axmy(mesh->Nlocal, 1.0, mesh->o_invLMM, o_Lf);
  }

  // Power iteration for the spectral radius of L (sets the stable explicit
  // step) and a Rayleigh quotient for its sign convention.
  void setupSmoothing(double tau)
  {
    auto mesh = nrs->meshV;
    const dlong N = mesh->Nlocal;
    std::vector<dfloat> v0(N);
    std::mt19937_64 gen(7 + platform->comm.mpiRank());
    std::uniform_real_distribution<double> uni(-1.0, 1.0);
    for (auto &v : v0) {
      v = uni(gen);
    }
    auto o_v = platform->device.malloc<dfloat>(N);
    auto o_w = platform->device.malloc<dfloat>(N);
    o_v.copyFrom(v0.data(), N);
    double lambda = 0, rayleigh = 0;
    for (int it = 0; it < 60; it++) {
      const double nv = std::sqrt(platform->linAlg->innerProd(N, o_v, o_v, comm()));
      platform->linAlg->scale(N, 1.0 / nv, o_v);
      applyLaplacian(o_v, o_w);
      rayleigh = platform->linAlg->innerProd(N, o_v, o_w, comm());
      lambda = std::sqrt(platform->linAlg->innerProd(N, o_w, o_w, comm()));
      o_v.copyFrom(o_w, N);
    }
    lapSign = rayleigh < 0 ? 1 : -1; // diffusion is f += dt * (-K) f
    dtSmooth = 1.5 / lambda;
    nSmooth = std::max(1, static_cast<int>(std::ceil(tau / dtSmooth)));
    dtSmooth = tau / nSmooth;
    if (platform->comm.mpiRank() == 0) {
      printf("elgpu: |lambda_max(L)| %.3e, Rayleigh %+.2e, filter as %d explicit steps of %.3e per field\n",
             lambda, rayleigh, nSmooth, dtSmooth);
    }
  }

  void diffuse(occa::memory &o_f)
  {
    auto o_Lf = platform->deviceMemoryPool.reserve<dfloat>(nrs->meshV->Nlocal);
    applyLaplacian(o_f, o_Lf);
    platform->linAlg->axpby(nrs->meshV->Nlocal, lapSign * dtSmooth, o_Lf, 1.0, o_f);
  }

  void updateSources(const double *dt)
  {
    auto mesh = nrs->meshV;
    auto &fluid = nrs->fluid;
    const auto offset = fluid->fieldOffset;
    const int n = particles->numParticles();
    const int poff = particles->fieldOffset();

    if (!o_q.isInitialized() || o_q.length() < static_cast<size_t>(nAcc * std::max(poff, 1))) {
      o_q = platform->device.malloc<dfloat>(nAcc * std::max(poff, 1));
      o_el = platform->device.malloc<dlong>(std::max(poff, 1));
      o_r = platform->device.malloc<dfloat>(3 * std::max(poff, 1));
      o_Up = platform->device.malloc<dfloat>(3 * std::max(poff, 1));
    }
    platform->linAlg->fill(nAcc * offset, 0.0, o_acc);
    // the current fluid velocity at the new positions (lpm_t's own
    // interpolate() would use the state extrapolated to the last RK stage);
    // collective, so every rank calls it, particles or not
    particles->interpolator().eval(3, offset, fluid->o_U, poff, o_Up);
    if (n) {
      forceKernel(n, poff, particles->dofId("velocity"), cfg.mu, cfg.rhof, particles->getDOF("x"), particles->getProp("dp"), particles->getProp("np"), o_Up, o_q);

      // after migrate() every particle is local: findpts' element and
      // reference coordinates are in particle order on the host
      auto &data = particles->interpolator().data();
      o_el.copyFrom(data.el.data(), n);
      o_r.copyFrom(data.r.data(), 3 * n);
      depositKernel(n, poff, mesh->Nq, offset, o_el, o_r, mesh->o_gllz, nAcc, o_q, o_acc);
    }
    oogs::startFinish(o_acc, nAcc, offset, ogsDfloat, ogsAdd, mesh->oogs);
    const bool dbg = getenv("ELGPU_DEBUG") != nullptr;
    for (int c = 0; c < nAcc; c++) {
      auto o_f = o_acc.slice(c * offset, mesh->Nlocal);
      platform->linAlg->axmy(mesh->Nlocal, 1.0, mesh->o_invLMM, o_f);
      double e0 = 0, e1 = 0;
      if (dbg && c == 3) {
        e0 = platform->linAlg->innerProd(mesh->Nlocal, o_f, o_f, comm());
      }
      for (int s = 0; s < nSmooth; s++) {
        diffuse(o_f);
      }
      if (dbg && c == 3) {
        e1 = platform->linAlg->innerProd(mesh->Nlocal, o_f, o_f, comm());
        if (platform->comm.mpiRank() == 0) {
          printf("elgpu dbg: phi l2^2 before/after smoothing %.3e %.3e\n", e0, e1);
        }
      }
    }

    auto o_phif = platform->deviceMemoryPool.reserve<dfloat>(mesh->Nlocal);
    platform->linAlg->axpbyz(mesh->Nlocal, -1.0, o_acc.slice(3 * offset, mesh->Nlocal), 1.0, o_one, o_phif);
    auto o_grad = opSEM::strongGrad(mesh, offset, o_phif);
    auto o_divFlux = opSEM::strongDivergence(mesh, offset, o_acc.slice(4 * offset, 3 * offset));
    sourceKernel(mesh->Nlocal, offset, cfg.qtlFromFlux ? 1 : 0, dt ? 1.0 / *dt : 0.0, o_acc, o_phifLast, o_divFlux,
                 fluid->o_U, o_grad, o_FU, o_div);
    if (!cfg.twoWay) {
      platform->linAlg->fill(fluid->fieldOffsetSum, 0.0, o_FU);
    }
    if (!dt || !cfg.qtl || !cfg.twoWay) {
      platform->linAlg->fill(offset, 0.0, o_div);
    } else if (cfg.zeroMeanDiv) {
      const double mean = platform->linAlg->innerProd(mesh->Nlocal, mesh->o_LMM, o_div, comm()) / volume;
      platform->linAlg->add(mesh->Nlocal, -mean, o_div);
    }
    o_phifLast.copyFrom(o_phif, mesh->Nlocal);

    if (getenv("ELGPU_DEBUG")) {
      // every rank takes part in each reduction, particles or not
      auto amax = [&](const occa::memory &o_a, dlong N) {
        double m = 0;
        if (N > 0) {
          std::vector<dfloat> a(N);
          o_a.copyTo(a.data(), N);
          for (auto v : a) {
            m = std::max(m, std::abs(static_cast<double>(v)));
          }
        }
        MPI_Allreduce(MPI_IN_PLACE, &m, 1, MPI_DOUBLE, MPI_MAX, comm());
        return m;
      };
      const double mq = amax(o_q, n ? 3 * poff : 0);
      const double mv = amax(particles->getDOF("velocity"), n ? 3 * poff : 0);
      const double mu = amax(o_Up, n ? 3 * poff : 0);
      double macc[3], mFU[3];
      for (int k = 0; k < 3; k++) {
        macc[k] = amax(o_acc.slice(k * offset, mesh->Nlocal), mesh->Nlocal);
        mFU[k] = amax(o_FU.slice(k * offset, mesh->Nlocal), mesh->Nlocal);
      }
      const double mphi = amax(o_acc.slice(3 * offset, mesh->Nlocal), mesh->Nlocal);
      const double mdiv = amax(o_div, mesh->Nlocal);
      if (platform->comm.mpiRank() == 0) {
        printf("elgpu dbg: n %d max|F_p| %.3e max|v_p| %.3e max|u_p| %.3e | max|F| %.2e %.2e %.2e max|FU| %.2e %.2e %.2e max phi %.3e max|div| %.3e\n",
               n, mq, mv, mu, macc[0], macc[1], macc[2], mFU[0], mFU[1], mFU[2], mphi, mdiv);
      }
    }
  }
};

} // namespace elgpu
