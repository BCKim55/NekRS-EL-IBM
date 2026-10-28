#!/usr/bin/env bash
# copy (not symlink) the UDF and .usr of the base case; the mesh and par are here
cd "$(dirname "$0")" && cp ../launch.udf ../launch.usr .
