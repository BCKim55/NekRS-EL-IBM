#!/usr/bin/env bash
# copy (not symlink) the mesh and .usr of the base case; the par and UDF are here
cd "$(dirname "$0")" && cp ../launch.{box,re2,usr} .
