#!/usr/bin/env bash
# copy (not symlink) the Re 100 hires mesh; everything else is in this directory
cd "$(dirname "$0")" && cp ../ibm_sphere_re100/hires/sphere.{box,re2} .
