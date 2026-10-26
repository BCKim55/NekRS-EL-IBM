#!/usr/bin/env bash
# copy (not symlink: nekRS's mesh repair diverges on symlinked inputs) the
# shared files of the base case, then run as usual
cd "$(dirname "$0")" && cp ../ibm_fixed_sphere/sphere.{box,re2,usr,udf} .
