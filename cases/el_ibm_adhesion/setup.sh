#!/usr/bin/env bash
# copy (not symlink) the shared files of el_ibm_collision; this directory holds the par
cd "$(dirname "$0")" && cp ../el_ibm_collision/collision.{box,re2,udf} .
