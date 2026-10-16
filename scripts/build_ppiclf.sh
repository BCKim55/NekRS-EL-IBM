#!/usr/bin/env bash
# Build a case-specific libppiclF.a. ppiclF compiles the user RHS and the
# array sizes (PPICLF_USER.h) into the library, so every case gets its own.
#
#   scripts/build_ppiclf.sh <case user_routines dir> <build dir>
#
# Sources come from the pinned external/ppiclF plus patches/ppiclF/*.patch.
# Writes <build dir>/udf.env with NEKRS_UDF_INCLUDES / NEKRS_UDF_LIBS for the UDF.
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
user=$(cd "$1" && pwd)
mkdir -p "$2"
out=$(cd "$2" && pwd)

FC=${FC:-mpif77}
CC=${CC:-mpicc}

src="$out/ppiclF"
stamp=$(cat "$user/PPICLF_USER.h" "$user/ppiclf_user.f" "$root/src/ppiclf_cbind.f" \
          "$root"/patches/ppiclF/*.patch | shasum | cut -c1-16)

if [ ! -f "$src/source/libppiclF.a" ] || [ "$(cat "$src/.stamp" 2>/dev/null)" != "$stamp" ]; then
  # fresh upstream sources; keep the downloaded/built gslib
  rsync -a --delete --exclude .git --exclude /.stamp \
    --exclude /3rd_party/gslib/gslib/ --exclude /3rd_party/gslib/lib/ \
    --exclude /3rd_party/gslib/include/ --exclude '/3rd_party/gslib/*.tar.gz' \
    "$root/external/ppiclF/" "$src/"
  for p in "$root"/patches/ppiclF/*.patch; do
    patch -s -p1 -d "$src" < "$p"
  done
  cp "$user/PPICLF_USER.h" "$src/source/"
  cat "$user/ppiclf_user.f" "$root/src/ppiclf_cbind.f" > "$src/source/ppiclf_user.f"

  # -fPIC: the library ends up inside the nekRS UDF shared object.
  CFLAGS="-O2 -fPIC" make -C "$src" FC="$FC" CC="$CC" FFLAGS="-cpp -std=legacy -O2 -fPIC"
  echo "$stamp" > "$src/.stamp"
fi

fortran_libs=$("$FC" --showme:link 2>/dev/null | tr ' ' '\n' | grep -E '^-(L|l)' | paste -sd ';' - || true)
gfortran_dir=$(dirname "$("$FC" -print-file-name=libgfortran.a)")
cat > "$out/udf.env" <<EOF
export NEKRS_UDF_INCLUDES="$root/src;$user;$src/source"
export NEKRS_UDF_LIBS="$src/source/libppiclF.a;$fortran_libs;-L$gfortran_dir;-lgfortran"
EOF
echo "ppiclF built: $src/source/libppiclF.a"
