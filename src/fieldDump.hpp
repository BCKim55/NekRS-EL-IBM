// Gather the GLL points of one coordinate plane (|coord - value| < tol) with
// velocity and pressure to a CSV on rank 0, for contour plots of a case's
// final state (scripts/plot_slices.py).
#pragma once
#include <mpi.h>
#include <cstdio>
#include <string>
#include <vector>

inline void dumpPlane(nrs_t *nrs, const std::string &name, int axis, double value, double tol)
{
  auto mesh = nrs->meshV;
  const dlong N = mesh->Nlocal;
  const auto off = nrs->fluid->fieldOffset;
  std::vector<dfloat> x(N), y(N), z(N), U(3 * off), P(off);
  mesh->o_x.copyTo(x.data(), N);
  mesh->o_y.copyTo(y.data(), N);
  mesh->o_z.copyTo(z.data(), N);
  nrs->fluid->o_U.copyTo(U.data(), 3 * off);
  nrs->fluid->o_P.copyTo(P.data(), off);
  const dfloat *c[3] = {x.data(), y.data(), z.data()};
  std::vector<double> rows;
  for (dlong i = 0; i < N; i++) {
    if (std::fabs(c[axis][i] - value) < tol) {
      rows.insert(rows.end(), {x[i], y[i], z[i], U[i], U[off + i], U[2 * off + i], P[i]});
    }
  }
  auto comm = platform->comm.mpiComm();
  const int size = platform->comm.mpiCommSize();
  int n = rows.size();
  std::vector<int> counts(size), displs(size, 0);
  MPI_Gather(&n, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);
  for (int r = 1; r < size; r++) {
    displs[r] = displs[r - 1] + counts[r - 1];
  }
  const int total = displs.back() + counts.back();
  std::vector<double> all(platform->comm.mpiRank() == 0 ? total : 0);
  MPI_Gatherv(rows.data(), n, MPI_DOUBLE, all.data(), counts.data(), displs.data(), MPI_DOUBLE, 0, comm);
  if (platform->comm.mpiRank() == 0) {
    FILE *f = fopen(name.c_str(), "w");
    fprintf(f, "x,y,z,u,v,w,p\n");
    for (int i = 0; i < total; i += 7) {
      fprintf(f, "%.6e,%.6e,%.6e,%.8e,%.8e,%.8e,%.8e\n", all[i], all[i + 1], all[i + 2], all[i + 3], all[i + 4], all[i + 5], all[i + 6]);
    }
    fclose(f);
  }
}
