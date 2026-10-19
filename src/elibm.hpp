#pragma once

// Coupled EL point particles and IBM rigid spheres for a nekRS UDF.
//
// ppiclF (src/ppiclf_user_elibm.f) carries two kinds of particles: EL point
// particles and the surface markers of every IBM body. The bodies themselves
// are few, so each rank keeps all of them here and integrates them
// identically. One call to step() after each fluid step:
//
//   1. particles: ppiclF RK3 with the new fluid velocity, optionally in
//      `substeps` sub-steps (contact stiffness usually needs a smaller dt
//      than the fluid); markers move with their bodies and compute the
//      direct-forcing force f = rho_f dV (u_body - u)/dt (as NekIBM);
//   2. bodies: (rho_b - rho_f) V dU/dt = -sum f_marker + F_contact
//      + (rho_b - rho_f) V g, likewise for the angular velocity with
//      -sum r x f_marker (Uhlmann 2005: the fluid inside the body is taken
//      to move rigidly; explicit, so rho_b/rho_f must exceed
//      1 + 3 sqrt(2 pi) sigma / R, checked in initialize()); markers are
//      then snapped to the new rigid positions;
//   3. projection: FU = F_EL / (1 - phi_p) + F_IBM goes to
//      fluid->o_explicitTerms() (userSource) and, with qtl, the volume-filtered
//      continuity div u = -(1/phi_f) D phi_f/Dt goes to fluid->o_div
//      (userDivergence) for the next fluid step.
//
// nekRS extrapolates explicit terms in time (EXT2: 2 E^n - E^{n-1}) to
// t^{n+1}, which is right for the smooth EL drag. For the IBM direct forcing
// it doubles the controller gain and makes it unstable, so applySource()
// writes the IBM part as E^n = (F - sum_{s>=1} c_s E^{n-s}) / c_0 and the
// step applies exactly F (marker slip then decays as (1 - gain)^n).
//
// Host-staged and CPU-bound (ppiclF); see docs/design.md for the GPU plan.

#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include "opSEM.hpp"
#include "ppiclfCoupler.hpp"

extern "C" {
void ppiclc_elibm_SetFluid(double *mu, double *rho, double *grav, double *dt);
void ppiclc_elibm_SetContact(double *kn, double *restitution, int *ifnb);
void ppiclc_elibm_SetBodies(int *nb,
                            double *x,
                            double *v,
                            double *w,
                            double *rot,
                            double *radius,
                            double *mass);
void ppiclc_elibm_PlaceMarkers();
void ppiclc_elibm_BeginStep();
}

namespace elibm
{

using vec3 = std::array<double, 3>;

struct Body {
  double radius = 0;
  double rho = 0;
  vec3 x{};
  vec3 v{};
  vec3 w{};
  std::array<double, 9> rot{1, 0, 0, 0, 1, 0, 0, 0, 1}; // body -> lab, column-major
  double volume() const { return 4.0 / 3.0 * M_PI * radius * radius * radius; }
};

struct ElParticle {
  vec3 x{};
  vec3 v{};
  double dp = 0;
  double rho = 0;
};

struct Config {
  double mu = 0;   // fluid dynamic viscosity
  double rhof = 0; // fluid density
  vec3 grav{};
  std::array<bool, 3> periodic{};
  double filterWidth = 0;   // Gaussian FWHM for projection and spreading
  double markerSpacing = 0; // ~ grid spacing; sets the marker count
  double kn = 0;            // contact stiffness (0: no contact)
  double restitution = 0.9;
  double neighborWidth = 0; // > 0 enables EL-EL contact (>= largest dp)
  int substeps = 1;
  int ioStep = 0; // ppiclF VTU interval in fluid steps (0: none)
  bool qtl = true;
  bool zeroMeanDiv = true; // closed or fully periodic domain
};

// Totals over all ranks, for diagnostics.
struct Totals {
  vec3 fluidMomentum{}; // rho_f int u dV (includes the fluid inside bodies)
  vec3 elMomentum{};
  vec3 bodyMomentum{}; // (rho_b - rho_f) V U, see header comment
  vec3 elForceOnFluid{};
  vec3 elForceProjected{}; // int F_EL dV
  vec3 ibmForceOnFluid{};
  vec3 ibmForceProjected{}; // int F_IBM dV
  vec3 source{};            // int FU dV, the momentum source of the next fluid step
  double divRms = 0;        // RMS of the imposed divergence (qtl)
  double divErrRms = 0;     // RMS of div u minus the imposed divergence
  double elMass = 0;
  long long nEl = 0;
  long long nMarker = 0;
};

class Coupling
{
public:
  Coupling(nrs_t *nrs, const Config &cfg) : nrs(nrs), cfg(cfg)
  {
    auto mesh = nrs->meshV;
    massLocal.resize(mesh->Nlocal);
    mesh->o_LMM.copyTo(massLocal.data(), mesh->Nlocal);

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

    o_FU.resize(nrs->fluid->fieldOffsetSum);
    o_div.resize(nrs->fluid->fieldOffset);
    nrs->userSource = [this](double) { applySource(); };
    if (cfg.qtl) {
      nrs->userDivergence = [this](double) { this->nrs->fluid->o_div.copyFrom(o_div, this->nrs->meshV->Nlocal); };
    }
  }

  const vec3 &domainLo() const { return lo; }
  const vec3 &domainHi() const { return hi; }
  std::vector<Body> &bodies() { return bodyList; }
  double slipRms() const { return slipRmsLast; }

  // Particles given on any rank are created on that rank; ppiclF bins them.
  void addEl(const ElParticle &p) { elList.push_back(p); }
  void addBody(const Body &b) { bodyList.push_back(b); }

  void initialize()
  {
    auto mesh = nrs->meshV;
    // Explicit body update: the fluid shell the markers drive (area x
    // Gaussian thickness) must be lighter than the body's excess mass, or a
    // body velocity change is fed back amplified every step. For a sphere,
    // rho_b / rho_f > 1 + 3 sqrt(2 pi) sigma / R.
    const double sigma = cfg.filterWidth / (2 * std::sqrt(2 * std::log(2.0)));
    for (const auto &b : bodyList) {
      const double minRatio = 1 + 3 * std::sqrt(2 * M_PI) * sigma / b.radius;
      nekrsCheck(b.rho <= minRatio * cfg.rhof,
                 comm(),
                 EXIT_FAILURE,
                 "IBM body density ratio %g must exceed %g for radius %g and filter width %g\n",
                 b.rho / cfg.rhof,
                 minRatio,
                 b.radius,
                 cfg.filterWidth);
    }

    ppiclf::initMPI(comm());

    std::vector<double> y0, rprop0;
    auto push = [&](const vec3 &x, const vec3 &v, std::array<double, PPICLF_LRP> &r) {
      y0.insert(y0.end(), {x[0], x[1], x[2], v[0], v[1], v[2]});
      rprop0.insert(rprop0.end(), r.begin(), r.end());
    };

    long long elOffset = elList.size();
    MPI_Exscan(MPI_IN_PLACE, &elOffset, 1, MPI_LONG_LONG, MPI_SUM, comm());
    if (platform->comm.mpiRank() == 0) {
      elOffset = 0;
    }
    for (size_t i = 0; i < elList.size(); i++) {
      const auto &p = elList[i];
      std::array<double, PPICLF_LRP> r{};
      r[PPICLF_R_JTYPE - 1] = PPICLF_TYPE_EL;
      r[PPICLF_R_JRHOP - 1] = p.rho;
      r[PPICLF_R_JDP - 1] = p.dp;
      r[PPICLF_R_JVOLP - 1] = M_PI / 6.0 * p.dp * p.dp * p.dp;
      r[PPICLF_R_JID - 1] = elOffset + i;
      push(p.x, p.v, r);
    }
    elList.clear();

    // markers of every body, created on rank 0
    if (platform->comm.mpiRank() == 0) {
      for (size_t b = 0; b < bodyList.size(); b++) {
        const auto &body = bodyList[b];
        const double area = 4 * M_PI * body.radius * body.radius;
        const int nm = std::max(12, static_cast<int>(std::lround(area / (cfg.markerSpacing * cfg.markerSpacing))));
        const double golden = M_PI * (3.0 - std::sqrt(5.0));
        // Marker volume = surface share x the Gaussian's normal thickness
        // sqrt(2 pi) sigma: a Gaussian sheet of that thickness integrates to
        // one across the surface, so spreading then interpolating returns
        // the marker's own velocity correction (gain 1). With the grid
        // spacing h instead (Uhlmann's choice for his narrower delta) the
        // gain would be h / (sqrt(2 pi) sigma) < 1 and slip would linger.
        const double thickness = std::sqrt(2 * M_PI) * sigma;
        for (int m = 0; m < nm; m++) {
          const double zr = 1.0 - (m + 0.5) * 2.0 / nm;
          const double rr = std::sqrt(1.0 - zr * zr);
          const vec3 ref = {body.radius * rr * std::cos(m * golden),
                            body.radius * rr * std::sin(m * golden),
                            body.radius * zr};
          std::array<double, PPICLF_LRP> r{};
          r[PPICLF_R_JTYPE - 1] = PPICLF_TYPE_MARKER;
          r[PPICLF_R_JVOLP - 1] = area / nm * thickness;
          r[PPICLF_R_JID - 1] = b;
          r[PPICLF_R_JRX - 1] = ref[0];
          r[PPICLF_R_JRY - 1] = ref[1];
          r[PPICLF_R_JRZ - 1] = ref[2];
          push(add(body.x, ref), body.v, r); // placed exactly below
        }
      }
    }

    int imethod = 1; // RK3
    int ndim = mesh->dim;
    int iendian = 0;
    int npart = y0.size() / PPICLF_LRS;
    y0.push_back(0); // never pass a null pointer
    rprop0.push_back(0);
    ppiclc_solve_InitParticle(&imethod, &ndim, &iendian, &npart, y0.data(), rprop0.data());

    double alpha = 1e-4; // Gaussian truncated where it falls to alpha of its peak
    int iwallm = 0;
    double filt = cfg.filterWidth;
    ppiclc_solve_InitGaussianFilter(&filt, &alpha, &iwallm);
    if (cfg.neighborWidth > 0) {
      double width = cfg.neighborWidth;
      ppiclc_solve_InitNeighborBin(&width);
    }
    ppiclf::initOverlapMesh(mesh);
    if (cfg.periodic[0]) {
      ppiclc_solve_InitPeriodicX(&lo[0], &hi[0]);
    }
    if (cfg.periodic[1]) {
      ppiclc_solve_InitPeriodicY(&lo[1], &hi[1]);
    }
    if (cfg.periodic[2]) {
      ppiclc_solve_InitPeriodicZ(&lo[2], &hi[2]);
    }

    double mu = cfg.mu, rhof = cfg.rhof, dt = nrs->dt[0];
    vec3 g = cfg.grav;
    ppiclc_elibm_SetFluid(&mu, &rhof, g.data(), &dt);
    double kn = cfg.kn, e = cfg.restitution;
    int ifnb = cfg.neighborWidth > 0 ? 1 : 0;
    ppiclc_elibm_SetContact(&kn, &e, &ifnb);
    pushBodies();
    ppiclc_elibm_PlaceMarkers();

    contactOnBody.assign(bodyList.size(), vec3{});
    updateSources(nullptr);
  }

  void step(double time, int tstep)
  {
    if (tstep == 0) {
      return;
    }
    auto mesh = nrs->meshV;
    auto &fluid = nrs->fluid;
    const double dt = nrs->dt[0];

    double mu = cfg.mu, rhof = cfg.rhof;
    vec3 g = cfg.grav;
    double dtf = dt;
    ppiclc_elibm_SetFluid(&mu, &rhof, g.data(), &dtf);

    // 1. particles
    ppiclc_elibm_BeginStep();
    double dts = dt / cfg.substeps;
    for (int s = 0; s < cfg.substeps; s++) {
      ppiclf::interpField(mesh, PPICLF_R_JUX, fluid->o_U.slice(0 * fluid->fieldOffset, mesh->Nlocal));
      ppiclf::interpField(mesh, PPICLF_R_JUY, fluid->o_U.slice(1 * fluid->fieldOffset, mesh->Nlocal));
      ppiclf::interpField(mesh, PPICLF_R_JUZ, fluid->o_U.slice(2 * fluid->fieldOffset, mesh->Nlocal));

      // ppiclF writes when mod(istep, iostep) == 0
      const bool out = (s == cfg.substeps - 1) && cfg.ioStep > 0;
      int istep = out ? tstep : -1;
      int iostep = out ? cfg.ioStep : (1 << 30);
      double t = time - dt + (s + 1) * dts;
      ppiclc_solve_IntegrateParticle(&istep, &iostep, &dts, &t);
    }

    // contact impulse each EL particle received from the bodies over the
    // step (RK-weighted in Fortran); the reaction goes to the nearest body
    std::vector<vec3> impulse(bodyList.size(), vec3{});
    if (bodyList.size()) {
      ppiclf::particles(yLoc, rLoc);
      for (int i = 0; i < ppiclf::numParticles(); i++) {
        const double *r = &rLoc[i * PPICLF_LRP];
        if (std::lround(r[PPICLF_R_JTYPE - 1]) != PPICLF_TYPE_EL) {
          continue;
        }
        const vec3 imp = {r[PPICLF_R_JIX - 1], r[PPICLF_R_JIY - 1], r[PPICLF_R_JIZ - 1]};
        if (imp[0] == 0 && imp[1] == 0 && imp[2] == 0) {
          continue;
        }
        const vec3 xp = {yLoc[i * PPICLF_LRS], yLoc[i * PPICLF_LRS + 1], yLoc[i * PPICLF_LRS + 2]};
        size_t bmin = 0;
        double gmin = std::numeric_limits<double>::max();
        for (size_t b = 0; b < bodyList.size(); b++) {
          const double gap = norm(minImage(sub(xp, bodyList[b].x))) - bodyList[b].radius;
          if (gap < gmin) {
            gmin = gap;
            bmin = b;
          }
        }
        for (int k = 0; k < 3; k++) {
          impulse[bmin][k] -= imp[k];
        }
      }
    }

    // 2. bodies
    if (bodyList.size()) {
      const size_t nb = bodyList.size();
      std::vector<double> buf(9 * nb, 0.0); // force, torque, contact impulse
      double slip2 = 0, nslip = 0;
      ppiclf::particles(yLoc, rLoc);
      for (int i = 0; i < ppiclf::numParticles(); i++) {
        const double *r = &rLoc[i * PPICLF_LRP];
        if (std::lround(r[PPICLF_R_JTYPE - 1]) != PPICLF_TYPE_MARKER) {
          continue;
        }
        const size_t b = std::lround(r[PPICLF_R_JID - 1]);
        const vec3 f = {-r[PPICLF_R_JFX - 1], -r[PPICLF_R_JFY - 1], -r[PPICLF_R_JFZ - 1]};
        const vec3 arm =
            minImage(sub(vec3{yLoc[i * PPICLF_LRS], yLoc[i * PPICLF_LRS + 1], yLoc[i * PPICLF_LRS + 2]}, bodyList[b].x));
        const vec3 t = cross(arm, f);
        for (int k = 0; k < 3; k++) {
          buf[9 * b + k] += f[k];
          buf[9 * b + 3 + k] += t[k];
        }
        const double s = norm(f) * dt / (cfg.rhof * r[PPICLF_R_JVOLP - 1]);
        slip2 += s * s;
        nslip += 1;
      }
      for (size_t b = 0; b < nb; b++) {
        for (int k = 0; k < 3; k++) {
          buf[9 * b + 6 + k] = impulse[b][k];
        }
      }
      MPI_Allreduce(MPI_IN_PLACE, buf.data(), buf.size(), MPI_DOUBLE, MPI_SUM, comm());
      double s2[2] = {slip2, nslip};
      MPI_Allreduce(MPI_IN_PLACE, s2, 2, MPI_DOUBLE, MPI_SUM, comm());
      slipRmsLast = s2[1] > 0 ? std::sqrt(s2[0] / s2[1]) : 0;

      for (size_t b = 0; b < nb; b++) {
        auto &body = bodyList[b];
        const double mEff = (body.rho - cfg.rhof) * body.volume();
        const double iEff = 0.4 * mEff * body.radius * body.radius;
        for (int k = 0; k < 3; k++) {
          contactOnBody[b][k] = buf[9 * b + 6 + k] / dt;
          body.v[k] += dt * ((buf[9 * b + k] + contactOnBody[b][k]) / mEff + cfg.grav[k]);
          body.x[k] += dt * body.v[k];
          body.w[k] += dt * buf[9 * b + 3 + k] / iEff;
          if (cfg.periodic[k]) {
            const double len = hi[k] - lo[k];
            body.x[k] = lo[k] + std::fmod(std::fmod(body.x[k] - lo[k], len) + len, len);
          }
        }
        rotate(body.rot, body.w, dt);
      }
      pushBodies();
      ppiclc_elibm_PlaceMarkers();
    }

    // 3. projection, momentum source, divergence
    updateSources(&dt);
  }

  Totals totals()
  {
    auto mesh = nrs->meshV;
    auto &fluid = nrs->fluid;
    Totals tot;

    for (int k = 0; k < 3; k++) {
      tot.fluidMomentum[k] =
          cfg.rhof * platform->linAlg->innerProd(mesh->Nlocal,
                                                 mesh->o_LMM,
                                                 fluid->o_U.slice(k * fluid->fieldOffset, mesh->Nlocal),
                                                 comm());
    }

    ppiclf::particles(yLoc, rLoc);
    std::vector<double> buf(14, 0.0);
    for (int i = 0; i < ppiclf::numParticles(); i++) {
      const double *r = &rLoc[i * PPICLF_LRP];
      const double *y = &yLoc[i * PPICLF_LRS];
      if (std::lround(r[PPICLF_R_JTYPE - 1]) == PPICLF_TYPE_EL) {
        const double m = r[PPICLF_R_JRHOP - 1] * r[PPICLF_R_JVOLP - 1];
        for (int k = 0; k < 3; k++) {
          buf[k] += m * y[PPICLF_JVX - 1 + k];
          buf[3 + k] += r[PPICLF_R_JFX - 1 + k];
        }
        buf[9] += m;
        buf[10] += 1;
      } else {
        for (int k = 0; k < 3; k++) {
          buf[6 + k] += r[PPICLF_R_JFX - 1 + k];
        }
        buf[11] += 1;
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, buf.data(), buf.size(), MPI_DOUBLE, MPI_SUM, comm());
    for (int k = 0; k < 3; k++) {
      tot.elMomentum[k] = buf[k];
      tot.elForceOnFluid[k] = buf[3 + k];
      tot.ibmForceOnFluid[k] = buf[6 + k];
      tot.elForceProjected[k] = intProjected[k];
      tot.ibmForceProjected[k] = intProjected[3 + k];
      tot.source[k] = intSource[k];
      for (const auto &b : bodyList) {
        tot.bodyMomentum[k] += (b.rho - cfg.rhof) * b.volume() * b.v[k];
      }
    }
    if (cfg.qtl) {
      auto o_divU = opSEM::strongDivergence(mesh, fluid->fieldOffset, fluid->o_U);
      std::vector<dfloat> divU(mesh->Nlocal), divQ(mesh->Nlocal);
      o_divU.copyTo(divU.data(), mesh->Nlocal);
      fluid->o_div.copyTo(divQ.data(), mesh->Nlocal);
      double d[3] = {0, 0, 0};
      for (dlong n = 0; n < mesh->Nlocal; n++) {
        d[0] += massLocal[n] * divQ[n] * divQ[n];
        d[1] += massLocal[n] * (divU[n] - divQ[n]) * (divU[n] - divQ[n]);
        d[2] += massLocal[n];
      }
      MPI_Allreduce(MPI_IN_PLACE, d, 3, MPI_DOUBLE, MPI_SUM, comm());
      tot.divRms = std::sqrt(d[0] / d[2]);
      tot.divErrRms = std::sqrt(d[1] / d[2]);
    }
    tot.elMass = buf[9];
    tot.nEl = std::llround(buf[10]);
    tot.nMarker = std::llround(buf[11]);
    return tot;
  }

private:
  nrs_t *nrs;
  Config cfg;
  vec3 lo{}, hi{};
  std::vector<ElParticle> elList;
  std::vector<Body> bodyList;
  std::vector<vec3> contactOnBody;
  std::vector<double> massLocal;
  std::vector<double> phifLast;
  std::vector<double> yLoc, rLoc;
  std::array<double, 6> intProjected{};
  vec3 intSource{};
  double slipRmsLast = 0;
  std::vector<dfloat> FUel, FUib;           // EL and IBM parts of the next source
  std::vector<std::vector<dfloat>> extHist; // IBM part as written, newest first
  deviceMemory<dfloat> o_FU;
  deviceMemory<dfloat> o_div;

  static MPI_Comm comm() { return platform->comm.mpiComm(); }

  static vec3 add(const vec3 &a, const vec3 &b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
  static vec3 sub(const vec3 &a, const vec3 &b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
  static double norm(const vec3 &a) { return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); }
  static vec3 cross(const vec3 &a, const vec3 &b)
  {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
  }

  vec3 minImage(vec3 d) const
  {
    for (int k = 0; k < 3; k++) {
      if (cfg.periodic[k]) {
        const double len = hi[k] - lo[k];
        d[k] -= len * std::nearbyint(d[k] / len);
      }
    }
    return d;
  }

  // rot <- exp([w dt]x) rot (Rodrigues)
  static void rotate(std::array<double, 9> &rot, const vec3 &w, double dt)
  {
    const double th = norm(w) * dt;
    if (th == 0) {
      return;
    }
    const vec3 a = {w[0] / norm(w), w[1] / norm(w), w[2] / norm(w)};
    const double c = std::cos(th), s = std::sin(th), t = 1 - c;
    const double R[3][3] = {{t * a[0] * a[0] + c, t * a[0] * a[1] - s * a[2], t * a[0] * a[2] + s * a[1]},
                            {t * a[0] * a[1] + s * a[2], t * a[1] * a[1] + c, t * a[1] * a[2] - s * a[0]},
                            {t * a[0] * a[2] - s * a[1], t * a[1] * a[2] + s * a[0], t * a[2] * a[2] + c}};
    std::array<double, 9> out{};
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          out[i + 3 * j] += R[i][k] * rot[k + 3 * j];
        }
      }
    }
    rot = out;
  }

  void applySource()
  {
    auto &fluid = nrs->fluid;
    std::vector<dfloat> c(fluid->o_coeffEXT.size());
    fluid->o_coeffEXT.copyTo(c.data(), c.size());

    std::vector<dfloat> E(FUib);
    for (size_t s = 1; s < c.size() && s <= extHist.size(); s++) {
      for (size_t n = 0; n < E.size(); n++) {
        E[n] -= c[s] * extHist[s - 1][n];
      }
    }
    for (auto &e : E) {
      e /= c[0];
    }
    extHist.insert(extHist.begin(), E);
    extHist.resize(std::min(extHist.size(), c.size() - 1));

    for (size_t n = 0; n < E.size(); n++) {
      E[n] += FUel[n];
    }
    o_FU.copyFrom(E);
    fluid->o_explicitTerms().copyFrom(o_FU, fluid->fieldOffsetSum);
  }

  void pushBodies()
  {
    int nb = bodyList.size();
    std::vector<double> x, v, w, rot, rad, mass;
    for (const auto &b : bodyList) {
      x.insert(x.end(), b.x.begin(), b.x.end());
      v.insert(v.end(), b.v.begin(), b.v.end());
      w.insert(w.end(), b.w.begin(), b.w.end());
      rot.insert(rot.end(), b.rot.begin(), b.rot.end());
      rad.push_back(b.radius);
      mass.push_back(b.rho * b.volume());
    }
    x.push_back(0);
    v.push_back(0);
    w.push_back(0);
    rot.push_back(0);
    rad.push_back(0);
    mass.push_back(0);
    ppiclc_elibm_SetBodies(&nb, x.data(), v.data(), w.data(), rot.data(), rad.data(), mass.data());
  }

  // dt == nullptr: initial call, D phi_f/Dt has no history yet
  void updateSources(const double *dt)
  {
    auto mesh = nrs->meshV;
    auto &fluid = nrs->fluid;
    const auto offset = fluid->fieldOffset;
    const dlong N = mesh->Nlocal;

    ppiclc_solve_Project();

    std::vector<double> phi;
    ppiclf::projectedField(mesh, PPICLF_P_JPHIP, phi);

    FUel.assign(fluid->fieldOffsetSum, 0.0);
    FUib.assign(fluid->fieldOffsetSum, 0.0);
    std::vector<double> f;
    for (int c = 0; c < 6; c++) {
      ppiclf::projectedField(mesh, PPICLF_P_JFELX + c, f);
      double sum = 0;
      for (dlong n = 0; n < N; n++) {
        if (c < 3) {
          FUel[n + c * offset] = f[n] / (1.0 - phi[n]);
        } else {
          FUib[n + (c - 3) * offset] = f[n];
        }
        sum += massLocal[n] * f[n];
      }
      MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, comm());
      intProjected[c] = sum;
    }
    for (int k = 0; k < 3; k++) {
      double sum = 0;
      for (dlong n = 0; n < N; n++) {
        sum += massLocal[n] * (FUel[n + k * offset] + FUib[n + k * offset]);
      }
      MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, comm());
      intSource[k] = sum;
    }

    if (!cfg.qtl) {
      return;
    }

    std::vector<double> phif(N);
    for (dlong n = 0; n < N; n++) {
      phif[n] = 1.0 - phi[n];
    }
    std::vector<dfloat> div(offset, 0.0);
    if (dt && phifLast.size() == phif.size()) {
      deviceMemory<dfloat> o_phif(std::vector<dfloat>(phif.begin(), phif.end()));
      auto o_grad = opSEM::strongGrad(mesh, offset, o_phif);
      std::vector<dfloat> grad(3 * offset), u(3 * offset);
      o_grad.copyTo(grad.data(), 3 * offset);
      fluid->o_U.copyTo(u.data(), 3 * offset);

      double sum[2] = {0, 0};
      for (dlong n = 0; n < N; n++) {
        double adv = 0;
        for (int k = 0; k < 3; k++) {
          adv += u[n + k * offset] * grad[n + k * offset];
        }
        div[n] = -((phif[n] - phifLast[n]) / *dt + adv) / phif[n];
        sum[0] += massLocal[n] * div[n];
        sum[1] += massLocal[n];
      }
      if (cfg.zeroMeanDiv) {
        MPI_Allreduce(MPI_IN_PLACE, sum, 2, MPI_DOUBLE, MPI_SUM, comm());
        for (dlong n = 0; n < N; n++) {
          div[n] -= sum[0] / sum[1];
        }
      }
    }
    phifLast = phif;
    o_div.copyFrom(div);
  }
};

} // namespace elibm
