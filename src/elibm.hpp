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
void ppiclc_elibm_SetContact(double *kn, double *restitution, int *ifnb, int *iattach, int *iplanar);
void ppiclc_elibm_SetBridge(double *sigma, double *theta, double *vbfrac, int *ibridgeEl, int *ibridgeIbm);
void ppiclc_elibm_SetBodies(int *nb,
                            double *x,
                            double *v,
                            double *w,
                            double *rot,
                            double *radius,
                            double *mass,
                            int *shape);
void ppiclc_elibm_PlaceMarkers();
void ppiclc_elibm_BeginStep();
void ppiclc_elibm_GetBridgeImpulse(int *nb, double *bimp);
void ppiclc_solve_InitParticleFilter(int *jfilt); // ppiclF patch 0006
}

namespace elibm
{

using vec3 = std::array<double, 3>;

struct Body {
  double radius = 0;
  double rho = 0;
  // > 0: a cylinder with its axis along z spanning this length (the periodic
  // slab thickness for a quasi-2D case): markers on the lateral surface,
  // contact on the in-plane distance, motion in x, y and spin about z only
  double length = 0;
  bool cylinder() const { return length > 0; }
  // a solid hemisphere of this radius: dome on the body +z side, flat face
  // at z = -3R/8 (x is the centre of mass; the sphere centre sits 3R/8
  // behind it on the axis); full rigid-body rotation with the tensor
  // I_zz = 2/5 m R^2, I_xx = I_yy = 83/320 m R^2
  bool hemisphere = false;
  int shape() const { return cylinder() ? 1 : hemisphere ? 2 : 0; } // as passed to Fortran
  static constexpr double kHemiOffset = 3.0 / 8.0; // sphere centre to centre of mass, in radii
  // orient the body axis (body z) along dir in the lab
  void setAxis(vec3 dir)
  {
    const double n = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
    const vec3 z = {dir[0] / n, dir[1] / n, dir[2] / n};
    const vec3 seed = std::abs(z[0]) < 0.9 ? vec3{1, 0, 0} : vec3{0, 1, 0};
    vec3 x = {seed[1] * z[2] - seed[2] * z[1], seed[2] * z[0] - seed[0] * z[2], seed[0] * z[1] - seed[1] * z[0]};
    const double nx = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    x = {x[0] / nx, x[1] / nx, x[2] / nx};
    const vec3 y = {z[1] * x[2] - z[2] * x[1], z[2] * x[0] - z[0] * x[2], z[0] * x[1] - z[1] * x[0]};
    rot = {x[0], x[1], x[2], y[0], y[1], y[2], z[0], z[1], z[2]}; // columns: body axes in the lab
  }
  vec3 axis() const { return {rot[6], rot[7], rot[8]}; }
  bool fixed = false; // held in place: markers still force the fluid, the body does not move
  // driven: the UDF sets x, v (and w) for the new time before step(); the
  // markers enforce that motion and hydroForceOnBodies() still reports the
  // fluid force (sum of the marker reactions; for an accelerating body the
  // force on the body proper is that plus rho_f V dU/dt, the internal fluid)
  bool prescribed = false;
  vec3 x{};
  vec3 v{};
  vec3 w{};
  std::array<double, 9> rot{1, 0, 0, 0, 1, 0, 0, 0, 1}; // body -> lab, column-major
  double volume() const
  {
    return cylinder() ? M_PI * radius * radius * length : (hemisphere ? 2.0 / 3.0 : 4.0 / 3.0) * M_PI * radius * radius * radius;
  }
};

struct ElParticle {
  vec3 x{};
  vec3 v{};
  double dp = 0;
  double rho = 0;
  double weight = 1; // droplets per parcel
};

struct Config {
  double mu = 0;   // fluid dynamic viscosity
  double rhof = 0; // fluid density
  vec3 grav{};
  std::array<bool, 3> periodic{};
  double filterWidth = 0;   // Gaussian FWHM of the EL projection (force, phi_p)
  double phiMax = 0;        // > 0: cap the projected phi_p in 1/(1-phi_p) and the continuity source (0: no cap)
  double filterWidthIbm = 0; // Gaussian FWHM of the IBM marker spreading (0: filterWidth)
  double markerSpacing = 0; // ~ grid spacing; sets the marker count
  double kn = 0;            // contact stiffness (0: no contact)
  double restitution = 0.9;
  double neighborWidth = 0; // > 0 enables EL-EL contact (>= largest dp)
  // EL particles stick to the first body they touch (perfectly inelastic),
  // then ride on it, exert no drag and add their mass to the body
  bool attachOnContact = false;
  bool planar = false; // EL particles keep their z (quasi-2D slab): no z velocity, no z force on the fluid
  // pendular liquid-bridge attraction (Willett et al. 2000) between EL
  // particles (bridgeEl) and EL particles and bodies (bridgeIbm): surface
  // tension sigma, contact angle (rad), bridge volume as a fraction of the
  // smaller sphere's volume; the bridge breaks at the rupture gap
  // (1 + theta/2) V^(1/3). Exclusive with attachOnContact on the body side.
  bool bridgeEl = false, bridgeIbm = false;
  double sigma = 0, contactAngle = 0, bridgeVolumeFraction = 0.1;
  double bridgeRupture(double dp) const // rupture gap of a bridge between spheres of diameter dp
  {
    return (1 + 0.5 * contactAngle) * std::cbrt(bridgeVolumeFraction * M_PI / 6 * dp * dp * dp);
  }
  int substeps = 1;
  bool interpEverySubstep = false; // re-interpolate the (frozen) fluid velocity at every substep, not only the first and last
  int ioStep = 0; // ppiclF VTU interval in fluid steps (0: none)
  bool twoWay = true; // false: particles feel the fluid, the fluid not the particles
  bool qtl = true;
  bool zeroMeanDiv = true; // closed or fully periodic domain
  // checkpoint slot for the projected particle volume fraction phi_p ("" : not written);
  // nekRS's field format only knows scalarNN, so pick a slot the case's own scalars do not use
  std::string phiCheckpointField = "scalar00";
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
  long long nAttached = 0;      // EL particles stuck to bodies
  double attachedMass = 0;      // their mass
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
    if (!cfg.phiCheckpointField.empty()) {
      // sized Nlocal, not fieldOffset: the checkpoint writer rejects padded fields
      o_phi.resize(mesh->Nlocal);
      platform->linAlg->fill(mesh->Nlocal, 0.0, o_phi);
      nrs->addUserCheckpointField(cfg.phiCheckpointField, {o_phi});
    }
    if (cfg.qtl) {
      nrs->userDivergence = [this](double) { this->nrs->fluid->o_div.copyFrom(o_div, this->nrs->meshV->Nlocal); };
    }
  }

  const vec3 &domainLo() const { return lo; }
  const vec3 &domainHi() const { return hi; }
  std::vector<Body> &bodies() { return bodyList; }
  double slipRms() const { return slipRmsLast; }
  // hydrodynamic force on each body (sum of the marker reactions) from the last step
  const std::vector<vec3> &hydroForceOnBodies() const { return hydroOnBody; }
  // contact force on each body from the EL parcels over the last step (spring-dashpot,
  // liquid bridges and the impulse of parcels that stuck), and its liquid-bridge part;
  // reported for fixed and prescribed bodies too
  const std::vector<vec3> &contactForceOnBodies() const { return contactOnBody; }
  const std::vector<vec3> &bridgeForceOnBodies() const { return bridgeOnBody; }
  // mass of the EL particles stuck to each body (attachOnContact)
  const std::vector<double> &attachedMass() const { return attachedMassOnBody; }

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
    const double ibmWidth = cfg.filterWidthIbm > 0 ? cfg.filterWidthIbm : cfg.filterWidth;
    const double sigma = ibmWidth / (2 * std::sqrt(2 * std::log(2.0)));
    for (const auto &b : bodyList) {
      const double minRatio = 1 + (b.cylinder() ? 2 : 3) * std::sqrt(2 * M_PI) * sigma / b.radius;
      nekrsCheck(!b.fixed && !b.prescribed && b.rho <= minRatio * cfg.rhof,
                 comm(),
                 EXIT_FAILURE,
                 "IBM body density ratio %g must exceed %g for radius %g and filter width %g\n",
                 b.rho / cfg.rhof,
                 minRatio,
                 b.radius,
                 ibmWidth);
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
      r[PPICLF_R_JVOLP - 1] = p.weight * M_PI / 6.0 * p.dp * p.dp * p.dp;
      r[PPICLF_R_JID - 1] = elOffset + i;
      r[PPICLF_R_JFILT - 1] = cfg.filterWidth;
      push(p.x, p.v, r);
    }
    elList.clear();
    {
      long long mine = y0.size() / PPICLF_LRS;
      MPI_Allreduce(&mine, &nElGlobal, 1, MPI_LONG_LONG, MPI_SUM, comm());
    }

    // markers of every body, created on rank 0
    if (platform->comm.mpiRank() == 0) {
      for (size_t b = 0; b < bodyList.size(); b++) {
        const auto &body = bodyList[b];
        const bool cyl = body.cylinder(), hemi = body.hemisphere;
        const double R = body.radius, ds = cfg.markerSpacing;
        const double area = cyl ? 2 * M_PI * R * body.length : (hemi ? 3 : 4) * M_PI * R * R;
        // cylinder: nz rings of nt markers on the lateral surface (no end caps: the slab is periodic)
        const int nt = std::max(12, static_cast<int>(std::lround(2 * M_PI * R / ds)));
        const int nz = std::max(1, static_cast<int>(std::lround(body.length / ds)));
        // sphere: Fibonacci points; hemisphere: the upper half of a sphere's
        // Fibonacci points (nd) plus rings on the flat face (nf), each set with
        // its own area share
        const int nsph = std::max(12, static_cast<int>(std::lround(4 * M_PI * R * R / (ds * ds))));
        const int nd = hemi ? (nsph + 1) / 2 : nsph;
        const int nr = std::max(1, static_cast<int>(std::lround(R / ds)));
        std::vector<vec3> flat; // flat-face markers, body frame (relative to the sphere centre)
        if (hemi) {
          for (int k = 0; k < nr; k++) {
            const double rk = (k + 0.5) * R / nr;
            const int nk = std::max(1, static_cast<int>(std::lround(2 * M_PI * rk / ds)));
            for (int j = 0; j < nk; j++) {
              const double th = 2 * M_PI * (j + 0.5 * (k % 2)) / nk;
              flat.push_back({rk * std::cos(th), rk * std::sin(th), 0.0});
            }
          }
        }
        const int nf = flat.size();
        const int nm = cyl ? nt * nz : nd + nf;
        const double golden = M_PI * (3.0 - std::sqrt(5.0));
        const vec3 com = {0, 0, hemi ? -Body::kHemiOffset * R : 0}; // centre of mass from the sphere centre
        // Marker volume = surface share x the Gaussian's normal thickness
        // sqrt(2 pi) sigma: a Gaussian sheet of that thickness integrates to
        // one across the surface, so spreading then interpolating returns
        // the marker's own velocity correction (gain 1). With the grid
        // spacing h instead (Uhlmann's choice for his narrower delta) the
        // gain would be h / (sqrt(2 pi) sigma) < 1 and slip would linger.
        const double thickness = std::sqrt(2 * M_PI) * sigma;
        for (int m = 0; m < nm; m++) {
          vec3 ref;
          double share = area / nm;
          if (cyl) {
            const double th = 2 * M_PI * (m % nt + 0.5 * ((m / nt) % 2)) / nt; // rings staggered by half a spacing
            ref = {R * std::cos(th), R * std::sin(th), (m / nt + 0.5) * body.length / nz - 0.5 * body.length};
          } else if (m < nd) {
            const double zr = 1.0 - (m + 0.5) * 2.0 / nsph; // hemisphere: m < nd keeps zr > 0
            const double rr = std::sqrt(1.0 - zr * zr);
            ref = {R * rr * std::cos(m * golden) - com[0], R * rr * std::sin(m * golden) - com[1], R * zr - com[2]};
            share = hemi ? 2 * M_PI * R * R / nd : area / nm;
          } else {
            ref = {flat[m - nd][0] - com[0], flat[m - nd][1] - com[1], flat[m - nd][2] - com[2]};
            share = M_PI * R * R / nf;
          }
          std::array<double, PPICLF_LRP> r{};
          r[PPICLF_R_JTYPE - 1] = PPICLF_TYPE_MARKER;
          r[PPICLF_R_JVOLP - 1] = share * thickness;
          r[PPICLF_R_JID - 1] = b;
          r[PPICLF_R_JRX - 1] = ref[0];
          r[PPICLF_R_JRY - 1] = ref[1];
          r[PPICLF_R_JRZ - 1] = ref[2];
          r[PPICLF_R_JFILT - 1] = ibmWidth;
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
    // ppiclF has one filter; its width sets the ghost and truncation distances,
    // so it is the larger of the two, and each particle projects with its own
    // (PPICLF_R_JFILT, patch 0006)
    double filt = std::max(cfg.filterWidth, ibmWidth);
    ppiclc_solve_InitGaussianFilter(&filt, &alpha, &iwallm);
    int jfilt = PPICLF_R_JFILT;
    ppiclc_solve_InitParticleFilter(&jfilt);
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
    int iattach = cfg.attachOnContact ? 1 : 0;
    int iplanar = cfg.planar ? 1 : 0;
    ppiclc_elibm_SetContact(&kn, &e, &ifnb, &iattach, &iplanar);
    if (cfg.bridgeEl || cfg.bridgeIbm) {
      if (cfg.bridgeIbm && cfg.attachOnContact) {
        throw std::runtime_error("elibm: bridgeIbm replaces attachOnContact; set one of them");
      }
      if (cfg.bridgeEl && !ifnb) {
        throw std::runtime_error("elibm: bridgeEl needs neighborWidth >= dp + bridgeRupture(dp)");
      }
      double sigma = cfg.sigma, theta = cfg.contactAngle, vb = cfg.bridgeVolumeFraction;
      int iel = cfg.bridgeEl ? 1 : 0, iibm = cfg.bridgeIbm ? 1 : 0;
      ppiclc_elibm_SetBridge(&sigma, &theta, &vb, &iel, &iibm);
    }
    pushBodies();
    ppiclc_elibm_PlaceMarkers();

    contactOnBody.assign(bodyList.size(), vec3{});
    bridgeOnBody.assign(bodyList.size(), vec3{});
    attachedMassOnBody.assign(bodyList.size(), 0.0);
    updateSources(nullptr);
  }

  // Add EL parcels during the run (the calling rank creates them; ppiclF
  // bins them). Collective: every rank must call it, with any count.
  void inject(const std::vector<ElParticle> &list)
  {
    long long offset = list.size();
    MPI_Exscan(MPI_IN_PLACE, &offset, 1, MPI_LONG_LONG, MPI_SUM, comm());
    if (platform->comm.mpiRank() == 0) {
      offset = 0;
    }
    std::vector<double> y, rprop;
    for (size_t i = 0; i < list.size(); i++) {
      const auto &p = list[i];
      std::array<double, PPICLF_LRP> r{};
      r[PPICLF_R_JTYPE - 1] = PPICLF_TYPE_EL;
      r[PPICLF_R_JRHOP - 1] = p.rho;
      r[PPICLF_R_JDP - 1] = p.dp;
      r[PPICLF_R_JVOLP - 1] = p.weight * M_PI / 6.0 * p.dp * p.dp * p.dp;
      r[PPICLF_R_JID - 1] = nElGlobal + offset + i;
      r[PPICLF_R_JFILT - 1] = cfg.filterWidth;
      y.insert(y.end(), {p.x[0], p.x[1], p.x[2], p.v[0], p.v[1], p.v[2]});
      rprop.insert(rprop.end(), r.begin(), r.end());
    }
    long long total = list.size();
    MPI_Allreduce(MPI_IN_PLACE, &total, 1, MPI_LONG_LONG, MPI_SUM, comm());
    nElGlobal += total;
    int n = list.size();
    y.push_back(0);
    rprop.push_back(0);
    ppiclc_solve_AddParticles(&n, y.data(), rprop.data());
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

    // prescribed bodies: the UDF has set their state for this step
    for (const auto &b : bodyList) {
      if (b.prescribed) {
        pushBodies();
        break;
      }
    }

    // 1. particles
    ppiclc_elibm_BeginStep();
    double dts = dt / cfg.substeps;
    // the fluid is frozen over the step: interpolate it to the particles
    // on the first and the last substep only (the first RK stage of the
    // substep consumes the queued fields; the last one gives the markers
    // the fluid velocity at their final position, as before) and project
    // the particles back once, on the last stage of the last substep
    for (int s = 0; s < cfg.substeps; s++) {
      if (cfg.interpEverySubstep || s == 0 || s == cfg.substeps - 1) {
        ppiclf::interpField(mesh, PPICLF_R_JUX, fluid->o_U.slice(0 * fluid->fieldOffset, mesh->Nlocal));
        ppiclf::interpField(mesh, PPICLF_R_JUY, fluid->o_U.slice(1 * fluid->fieldOffset, mesh->Nlocal));
        ppiclf::interpField(mesh, PPICLF_R_JUZ, fluid->o_U.slice(2 * fluid->fieldOffset, mesh->Nlocal));
      }
      int proj = (s == cfg.substeps - 1) ? 1 : 0;
      ppiclc_solve_SetProjStep(&proj);

      // ppiclF writes when mod(istep, iostep) == 0
      const bool out = (s == cfg.substeps - 1) && cfg.ioStep > 0;
      int istep = out ? tstep : -1;
      int iostep = out ? cfg.ioStep : (1 << 30);
      double t = time - dt + (s + 1) * dts;
      ppiclc_solve_IntegrateParticle(&istep, &iostep, &dts, &t);
    }

    // contact impulse each EL particle received from the bodies over the
    // step (RK-weighted in Fortran); the reaction goes to the nearest body
    std::vector<vec3> impulse(bodyList.size(), vec3{}), impulseTorque(bodyList.size(), vec3{});
    std::vector<double> mAtt(bodyList.size(), 0.0); // mass stuck to each body
    if (bodyList.size()) {
      ppiclf::particles(yLoc, rLoc);
      for (int i = 0; i < ppiclf::numParticles(); i++) {
        const double *r = &rLoc[i * PPICLF_LRP];
        if (std::lround(r[PPICLF_R_JTYPE - 1]) != PPICLF_TYPE_EL) {
          continue;
        }
        const int att = std::lround(r[PPICLF_R_JATT - 1]);
        if (att > 0) {
          mAtt[att - 1] += r[PPICLF_R_JRHOP - 1] * r[PPICLF_R_JVOLP - 1];
        }
        const vec3 imp = {r[PPICLF_R_JIX - 1], r[PPICLF_R_JIY - 1], r[PPICLF_R_JIZ - 1]};
        if (imp[0] == 0 && imp[1] == 0 && imp[2] == 0) {
          continue;
        }
        const vec3 xp = {yLoc[i * PPICLF_LRS], yLoc[i * PPICLF_LRS + 1], yLoc[i * PPICLF_LRS + 2]};
        size_t bmin = 0;
        double gmin = std::numeric_limits<double>::max();
        for (size_t b = 0; b < bodyList.size(); b++) {
          const auto &bb = bodyList[b];
          vec3 d = minImage(sub(xp, bb.x));
          if (bb.cylinder()) {
            d[2] = 0;
          }
          if (bb.hemisphere) { // distance from the sphere centre, good enough to pick the body
            const vec3 ax = bb.axis();
            for (int k = 0; k < 3; k++) {
              d[k] += Body::kHemiOffset * bb.radius * ax[k];
            }
          }
          const double gap = norm(d) - bb.radius;
          if (gap < gmin) {
            gmin = gap;
            bmin = b;
          }
        }
        // reaction on the body: the impulse and, for a hemisphere, its moment
        // about the centre of mass (a sphere's or cylinder's contact normal
        // passes through the axis: no torque, and the end-of-step arm would
        // only add roundoff)
        const vec3 arm = minImage(sub(xp, bodyList[bmin].x));
        const vec3 tq = bodyList[bmin].hemisphere ? cross(arm, imp) : vec3{};
        for (int k = 0; k < 3; k++) {
          impulse[bmin][k] -= imp[k];
          impulseTorque[bmin][k] -= tq[k];
        }
      }
    }

    // 2. bodies
    if (bodyList.size()) {
      const size_t nb = bodyList.size();
      constexpr int nbuf = 13; // force, torque, contact impulse, attached mass, contact torque impulse
      std::vector<double> buf(nbuf * nb, 0.0);
      std::vector<double> bimp(3 * nb);
      int nbi = nb;
      ppiclc_elibm_GetBridgeImpulse(&nbi, bimp.data());
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
          buf[nbuf * b + k] += f[k];
          buf[nbuf * b + 3 + k] += t[k];
        }
        const double s = norm(f) * dt / (cfg.rhof * r[PPICLF_R_JVOLP - 1]);
        slip2 += s * s;
        nslip += 1;
      }
      for (size_t b = 0; b < nb; b++) {
        for (int k = 0; k < 3; k++) {
          buf[nbuf * b + 6 + k] = impulse[b][k];
          buf[nbuf * b + 10 + k] = impulseTorque[b][k];
        }
        buf[nbuf * b + 9] = mAtt[b];
        for (int k = 0; k < 3; k++) {
          bimp[3 * b + k] = -bimp[3 * b + k]; // reaction on the body
        }
      }
      MPI_Allreduce(MPI_IN_PLACE, buf.data(), buf.size(), MPI_DOUBLE, MPI_SUM, comm());
      // reduced apart from buf so the body sums keep their reduction order (bit-identical free-body runs)
      MPI_Allreduce(MPI_IN_PLACE, bimp.data(), bimp.size(), MPI_DOUBLE, MPI_SUM, comm());
      double s2[2] = {slip2, nslip};
      MPI_Allreduce(MPI_IN_PLACE, s2, 2, MPI_DOUBLE, MPI_SUM, comm());
      slipRmsLast = s2[1] > 0 ? std::sqrt(s2[0] / s2[1]) : 0;

      hydroOnBody.resize(nb);
      for (size_t b = 0; b < nb; b++) {
        auto &body = bodyList[b];
        for (int k = 0; k < 3; k++) {
          hydroOnBody[b][k] = buf[nbuf * b + k];
          contactOnBody[b][k] = buf[nbuf * b + 6 + k] / dt;
          bridgeOnBody[b][k] = bimp[3 * b + k] / dt;
        }
        if (body.fixed) {
          continue;
        }
        if (body.prescribed) {
          for (int k = 0; k < 3; k++) {
            if (cfg.periodic[k]) {
              const double len = hi[k] - lo[k];
              body.x[k] = lo[k] + std::fmod(std::fmod(body.x[k] - lo[k], len) + len, len);
            }
          }
          continue;
        }
        const double mEff = (body.rho - cfg.rhof) * body.volume();
        const double mTot = mEff + buf[nbuf * b + 9]; // plus the stuck particles
        const bool cyl = body.cylinder();
        const double iEff = (cyl ? 0.5 : 0.4) * mEff * body.radius * body.radius;
        attachedMassOnBody[b] = buf[nbuf * b + 9];
        vec3 torque; // hydrodynamic + contact, about the centre of mass
        for (int k = 0; k < 3; k++) {
          torque[k] = buf[nbuf * b + 3 + k] + buf[nbuf * b + 10 + k] / dt;
        }
        for (int k = 0; k < 3; k++) {
          if (cyl && k == 2) { // a cylinder stays in its plane and spins about z only
            body.v[2] = body.w[0] = body.w[1] = 0;
            body.w[2] += dt * torque[2] / iEff;
            continue;
          }
          body.v[k] += dt * ((buf[nbuf * b + k] + contactOnBody[b][k] + mTot * cfg.grav[k]) / mTot);
          body.x[k] += dt * body.v[k];
          if (!cyl && !body.hemisphere) {
            body.w[k] += dt * torque[k] / iEff;
          }
          if (cfg.periodic[k]) {
            const double len = hi[k] - lo[k];
            body.x[k] = lo[k] + std::fmod(std::fmod(body.x[k] - lo[k], len) + len, len);
          }
        }
        if (body.hemisphere) {
          // Euler's equations in the body frame with I = diag(Ixx, Ixx, Izz)
          const double Ixx = 83.0 / 320.0 * mEff * body.radius * body.radius;
          const double Izz = 0.4 * mEff * body.radius * body.radius;
          const double I[3] = {Ixx, Ixx, Izz};
          vec3 wb{}, tb{};
          for (int i = 0; i < 3; i++) {
            for (int k = 0; k < 3; k++) { // rot^T
              wb[i] += body.rot[k + 3 * i] * body.w[k];
              tb[i] += body.rot[k + 3 * i] * torque[k];
            }
          }
          const vec3 Iw = {I[0] * wb[0], I[1] * wb[1], I[2] * wb[2]};
          const vec3 gyro = cross(wb, Iw);
          for (int i = 0; i < 3; i++) {
            wb[i] += dt * (tb[i] - gyro[i]) / I[i];
          }
          for (int k = 0; k < 3; k++) {
            body.w[k] = 0;
            for (int i = 0; i < 3; i++) {
              body.w[k] += body.rot[k + 3 * i] * wb[i];
            }
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
        if (std::lround(r[PPICLF_R_JATT - 1]) > 0) {
          buf[12] += 1;
          buf[13] += m;
        }
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
    tot.nAttached = std::llround(buf[12]);
    tot.attachedMass = buf[13];
    return tot;
  }

private:
  nrs_t *nrs;
  Config cfg;
  vec3 lo{}, hi{};
  std::vector<ElParticle> elList;
  long long nElGlobal = 0;
  std::vector<Body> bodyList;
  std::vector<vec3> contactOnBody, bridgeOnBody, hydroOnBody;
  std::vector<double> attachedMassOnBody;
  std::vector<double> massLocal;
  std::vector<double> phifLast;
  std::vector<double> yLoc, rLoc;
  std::array<double, 6> intProjected{};
  vec3 intSource{};
  // extrema of the projected fields on the last updateSources call (all ranks)
  std::array<double, 3> srcMax{}; // max phi_p, max |F_EL/(1-phi_p)|, max |div source|
  double slipRmsLast = 0;
  std::vector<dfloat> FUel, FUib;           // EL and IBM parts of the next source
  std::vector<std::vector<dfloat>> extHist; // IBM part as written, newest first
  deviceMemory<dfloat> o_FU;
  deviceMemory<dfloat> o_div;
  deviceMemory<dfloat> o_phi; // projected phi_p for the checkpoint files

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
    std::vector<int> shape;
    for (const auto &b : bodyList) {
      shape.push_back(b.shape());
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
    shape.push_back(0);
    ppiclc_elibm_SetBodies(&nb, x.data(), v.data(), w.data(), rot.data(), rad.data(), mass.data(), shape.data());
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
    if (cfg.phiMax > 0) {
      for (dlong n = 0; n < N; n++) {
        phi[n] = std::min(phi[n], cfg.phiMax);
      }
    }
    if (o_phi.size()) {
      o_phi.copyFrom(std::vector<dfloat>(phi.begin(), phi.begin() + N));
    }

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
    srcMax = {0, 0, 0};
    for (dlong n = 0; n < N; n++) {
      srcMax[0] = std::max(srcMax[0], phi[n]);
      for (int k = 0; k < 3; k++) {
        srcMax[1] = std::max(srcMax[1], std::abs(FUel[n + k * offset]));
      }
    }
    if (!cfg.twoWay) {
      std::fill(FUel.begin(), FUel.end(), 0.0);
      std::fill(FUib.begin(), FUib.end(), 0.0);
    }
    for (int k = 0; k < 3; k++) {
      double sum = 0;
      for (dlong n = 0; n < N; n++) {
        sum += massLocal[n] * (FUel[n + k * offset] + FUib[n + k * offset]);
      }
      MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, comm());
      intSource[k] = sum;
    }

    if (!cfg.qtl || !cfg.twoWay) {
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
        srcMax[2] = std::max(srcMax[2], std::abs((double)div[n]));
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

public:
  // max phi_p, max |F_EL/(1-phi_p)|, max |continuity source| of the last projection (collective)
  std::array<double, 3> sourceExtrema() const
  {
    auto m = srcMax;
    MPI_Allreduce(MPI_IN_PLACE, m.data(), 3, MPI_DOUBLE, MPI_MAX, comm());
    return m;
  }
};

} // namespace elibm
