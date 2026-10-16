#pragma once

// Host-staged coupling between nekRS (device fields) and ppiclF (CPU library).
//
// ppiclF sees nekRS's own element partition as its "overlap mesh": every rank
// passes its local GLL coordinates once, then before each particle step the
// fields to interpolate, and after it reads back the projected fields.
// nekRS stores a field as Nelements blocks of Nq^3 points (i fastest), which is
// exactly ppiclF's (LEX,LEY,LEZ,NEE) layout, so no reordering is needed.
//
// Include this from a .udf; the case's PPICLF_USER.h must be on the include path.

#include <mpi.h>
#include <type_traits>
#include <vector>

#include "PPICLC.h"

extern "C" {
// src/ppiclf_cbind.f
void ppiclc_solve_InterpFieldUser(int *jp, double *fld);
void ppiclc_solve_GetProFld(int *m, double *fld);
void ppiclc_solve_GetNumParticles(int *npart);
void ppiclc_solve_GetParticles(double *y, double *rprop);
}

static_assert(std::is_same<dfloat, double>::value, "ppiclF is real*8 only; build nekRS with double dfloat");

namespace ppiclf
{

inline void initMPI(MPI_Comm comm)
{
  int fcomm = MPI_Comm_c2f(comm);
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  ppiclc_comm_InitMPI(&fcomm, &rank, &size);
}

inline void initOverlapMesh(mesh_t *mesh)
{
  nekrsCheck(mesh->Nq != PPICLF_LEX || mesh->Nq != PPICLF_LEY || mesh->Nq != PPICLF_LEZ,
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "PPICLF_LEX/LEY/LEZ (%d) must equal polynomialOrder+1 (%d)\n",
             PPICLF_LEX,
             mesh->Nq);
  nekrsCheck(mesh->Nelements > PPICLF_LEE,
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "PPICLF_LEE (%d) < local element count (%d)\n",
             PPICLF_LEE,
             mesh->Nelements);

  std::vector<double> x(mesh->Nlocal), y(mesh->Nlocal), z(mesh->Nlocal);
  mesh->o_x.copyTo(x.data(), mesh->Nlocal);
  mesh->o_y.copyTo(y.data(), mesh->Nlocal);
  mesh->o_z.copyTo(z.data(), mesh->Nlocal);

  int nel = mesh->Nelements;
  int nq = mesh->Nq;
  ppiclc_comm_InitOverlapMesh(&nel, &nq, &nq, &nq, x.data(), y.data(), z.data());
}

// Hand one scalar field (Nlocal entries of o_fld, device) to ppiclF for
// interpolation into particle property jp (1-based PPICLF_R_* index).
// ppiclF copies it into its own storage, so the host buffer is reused.
inline void interpField(mesh_t *mesh, int jp, const occa::memory &o_fld)
{
  static std::vector<double> buf;
  buf.resize(mesh->Nlocal);
  o_fld.copyTo(buf.data(), mesh->Nlocal);
  ppiclc_solve_InterpFieldUser(&jp, buf.data());
}

// Copy projected field m (1-based PPICLF_P_* index) into o_fld (device).
inline void projectedField(mesh_t *mesh, int m, occa::memory &o_fld)
{
  static std::vector<double> buf;
  buf.resize(mesh->Nlocal);
  ppiclc_solve_GetProFld(&m, buf.data());
  o_fld.copyFrom(buf.data(), mesh->Nlocal);
}

inline int numParticles()
{
  int n;
  ppiclc_solve_GetNumParticles(&n);
  return n;
}

// Rank-local particle state, (PPICLF_LRS, n) and (PPICLF_LRP, n) column-major.
inline void particles(std::vector<double> &y, std::vector<double> &rprop)
{
  const int n = numParticles();
  y.resize(static_cast<size_t>(PPICLF_LRS) * n);
  rprop.resize(static_cast<size_t>(PPICLF_LRP) * n);
  ppiclc_solve_GetParticles(y.data(), rprop.data());
}

} // namespace ppiclf
