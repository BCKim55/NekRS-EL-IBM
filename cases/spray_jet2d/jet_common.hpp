// Shared by jet.udf (ppiclF path) and ../spray_jet2d_gpu/jet.udf (device
// path): parameters, the deterministic parcel generator, the momentum
// budget and the field dump, so that the two runs see identical input.
#pragma once

#include <random>

struct JetParams {
  double uJet, wSlot, xNozzle, phiInj, dp, rhop, filterWidth;
  int parcelsPerStep, reportInterval, qtl;
  void read()
  {
    auto get = [](const char *k, auto &v) { platform->par->extract("casedata", k, v); };
    get("u_jet", uJet);
    get("w_slot", wSlot);
    get("x_nozzle", xNozzle);
    get("phi_inj", phiInj);
    get("dp", dp);
    get("rhop", rhop);
    get("parcels_per_step", parcelsPerStep);
    get("filter_width", filterWidth);
    get("qtl", qtl);
    get("report_interval", reportInterval);
  }
  double dropletVolume() const { return M_PI / 6.0 * dp * dp * dp; }
  // liquid volume entering per step, and droplets per parcel
  double liquidPerStep(double dt, double lz) const { return uJet * wSlot * lz * phiInj * dt; }
  double weight(double dt, double lz) const { return liquidPerStep(dt, lz) / (parcelsPerStep * dropletVolume()); }
};

// parcels for one step, column-major (3, n), all on rank 0; the positions
// are spread over the distance travelled during the step so that the
// stream is continuous
struct JetInjector {
  const JetParams &p;
  double lz, yc, dt;
  std::mt19937_64 gen{2026};
  std::uniform_real_distribution<double> uni{0.0, 1.0};

  void parcels(std::vector<double> &x, std::vector<double> &v, std::vector<double> &d, std::vector<double> &rho, std::vector<double> &w)
  {
    x.clear(), v.clear(), d.clear(), rho.clear(), w.clear();
    if (platform->comm.mpiRank() != 0) {
      return;
    }
    const double np = p.weight(dt, lz);
    for (int i = 0; i < p.parcelsPerStep; i++) {
      x.insert(x.end(), {p.xNozzle + uni(gen) * p.uJet * dt, yc + (uni(gen) - 0.5) * p.wSlot, uni(gen) * lz});
      v.insert(v.end(), {p.uJet, 0.0, 0.0});
      d.push_back(p.dp);
      rho.push_back(p.rhop);
      w.push_back(np);
    }
  }
};

// gather rows of `ncol` doubles from every rank and write them as a CSV
inline void gatherRows(const std::vector<double> &rows, int ncol, const std::string &name)
{
  auto comm = platform->comm.mpiComm();
  const int size = platform->comm.mpiCommSize();
  int n = rows.size();
  std::vector<int> counts(size), displs(size, 0);
  MPI_Gather(&n, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);
  for (int r = 1; r < size; r++) {
    displs[r] = displs[r - 1] + counts[r - 1];
  }
  std::vector<double> all(platform->comm.mpiRank() == 0 ? displs.back() + counts.back() : 0);
  MPI_Gatherv(rows.data(), n, MPI_DOUBLE, all.data(), counts.data(), displs.data(), MPI_DOUBLE, 0, comm);
  if (platform->comm.mpiRank() == 0) {
    FILE *f = fopen(name.c_str(), "w");
    fprintf(f, "x,y,z,u\n");
    for (size_t i = 0; i < all.size(); i += ncol) {
      for (int c = 0; c < ncol; c++) {
        fprintf(f, "%.8e%s", all[i + c], c + 1 < ncol ? "," : "\n");
      }
    }
    fclose(f);
  }
}

// gather x, y, z, u, v, w of every GLL point to rank 0 and write a CSV
inline void dumpField(nrs_t *nrs, const std::string &name)
{
  auto mesh = nrs->meshV;
  const dlong N = mesh->Nlocal;
  const auto off = nrs->fluid->fieldOffset;
  std::vector<dfloat> buf(6 * N), U(3 * off);
  mesh->o_x.copyTo(buf.data(), N);
  mesh->o_y.copyTo(buf.data() + N, N);
  mesh->o_z.copyTo(buf.data() + 2 * N, N);
  nrs->fluid->o_U.copyTo(U.data(), 3 * off);
  for (int k = 0; k < 3; k++) {
    std::copy(U.begin() + k * off, U.begin() + k * off + N, buf.begin() + (3 + k) * N);
  }
  auto comm = platform->comm.mpiComm();
  const int size = platform->comm.mpiCommSize();
  int n = N;
  std::vector<int> counts(size), displs(size, 0);
  MPI_Gather(&n, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);
  for (int r = 1; r < size; r++) {
    displs[r] = displs[r - 1] + counts[r - 1];
  }
  const int total = displs.back() + counts.back();
  std::vector<std::vector<dfloat>> cols(6, std::vector<dfloat>(platform->comm.mpiRank() == 0 ? total : 0));
  for (int c = 0; c < 6; c++) {
    MPI_Gatherv(buf.data() + c * N, N, MPI_DOUBLE, cols[c].data(), counts.data(), displs.data(), MPI_DOUBLE, 0, comm);
  }
  if (platform->comm.mpiRank() == 0) {
    FILE *f = fopen(name.c_str(), "w");
    fprintf(f, "x,y,z,u,v,w\n");
    for (int i = 0; i < total; i++) {
      fprintf(f, "%.6e,%.6e,%.6e,%.8e,%.8e,%.8e\n", cols[0][i], cols[1][i], cols[2][i], cols[3][i], cols[4][i], cols[5][i]);
    }
    fclose(f);
  }
}
