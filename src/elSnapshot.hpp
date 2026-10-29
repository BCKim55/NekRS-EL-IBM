// EL parcel snapshots for ParaView (scripts/snapshots_to_vtk.py): gather the
// EL parcels to rank 0 and write snapshots/el_XXXXX.csv (id, x, y, z, u, v,
// w, contact, dp, attached body or 0) and append the first body to
// snapshots/body.csv (step, t, x, y, z, u, v, w, radius).
#pragma once
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <vector>

#include "elibm.hpp"

namespace elibm {

constexpr int kElCols = 10;

// rank-local EL parcels as rows of (id, x, y, z, u, v, w, contact, dp, attached)
inline std::vector<double> elRows()
{
  std::vector<double> y, r, rows;
  ppiclf::particles(y, r);
  for (int i = 0; i < ppiclf::numParticles(); i++) {
    const double *ri = &r[i * PPICLF_LRP];
    if (std::lround(ri[PPICLF_R_JTYPE - 1]) != PPICLF_TYPE_EL) {
      continue;
    }
    const double *yi = &y[i * PPICLF_LRS];
    const bool touching = ri[PPICLF_R_JCX - 1] != 0 || ri[PPICLF_R_JCY - 1] != 0 || ri[PPICLF_R_JCZ - 1] != 0;
    rows.insert(rows.end(), {ri[PPICLF_R_JID - 1], yi[0], yi[1], yi[2], yi[3], yi[4], yi[5], touching ? 1.0 : 0.0,
                             ri[PPICLF_R_JDP - 1], ri[PPICLF_R_JATT - 1]});
  }
  return rows;
}

// all EL rows (kElCols each) on rank 0, empty elsewhere; collective
inline std::vector<double> gatherElRows()
{
  auto comm = platform->comm.mpiComm();
  const auto rows = elRows();
  int n = rows.size();
  std::vector<int> counts(platform->comm.mpiCommSize()), displs(counts.size(), 0);
  MPI_Gather(&n, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);
  std::vector<double> all;
  if (platform->comm.mpiRank() == 0) {
    for (size_t r = 1; r < counts.size(); r++) {
      displs[r] = displs[r - 1] + counts[r - 1];
    }
    all.resize(displs.back() + counts.back());
  }
  MPI_Gatherv(rows.data(), n, MPI_DOUBLE, all.data(), counts.data(), displs.data(), MPI_DOUBLE, 0, comm);
  return all;
}

inline void writeSnapshot(Coupling &coupling, double time, int tstep)
{
  const auto all = gatherElRows();
  if (platform->comm.mpiRank() != 0) {
    return;
  }

  std::filesystem::create_directories("snapshots");
  char name[64];
  snprintf(name, sizeof(name), "snapshots/el_%05d.csv", tstep);
  FILE *f = fopen(name, "w");
  fprintf(f, "id,x,y,z,u,v,w,contact,dp,attached\n");
  for (size_t i = 0; i < all.size(); i += kElCols) {
    fprintf(f, "%d,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%d,%.6e,%d\n", static_cast<int>(all[i]), all[i + 1], all[i + 2],
            all[i + 3], all[i + 4], all[i + 5], all[i + 6], static_cast<int>(all[i + 7]), all[i + 8],
            static_cast<int>(std::lround(all[i + 9])));
  }
  fclose(f);

  const auto &b = coupling.bodies()[0];
  const bool first = tstep == 0;
  f = fopen("snapshots/body.csv", first ? "w" : "a");
  if (first) {
    fprintf(f, "step,t,x,y,z,u,v,w,radius\n");
  }
  fprintf(f, "%d,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e\n", tstep, time, b.x[0], b.x[1], b.x[2], b.v[0], b.v[1],
          b.v[2], b.radius);
  fclose(f);
}

} // namespace elibm
