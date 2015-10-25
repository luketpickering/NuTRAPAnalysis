#!/bin/bash

#If it was sourced as . setup.sh then you can't scrub off the end... assume that
#we are in the correct directory.
if ! echo "${BASH_SOURCE}" | grep "/" --silent; then
  SETUPDIR=$(readlink -f $PWD)
else
  SETUPDIR=$(readlink -f ${BASH_SOURCE%/*})
fi

if ! [[ ":$PATH:" == *":${SETUPDIR}/NuTRAPAnalysis/bin:"* ]]; then
  export PATH=${SETUPDIR}/NuTRAPAnalysis/bin:$PATH
fi
if ! [[ ":$LD_LIBRARY_PATH:" == *":${SETUPDIR}/NuTRAPAnalysis/lib:"* ]]; then
  export LD_LIBRARY_PATH=${SETUPDIR}/NuTRAPAnalysis/lib:$LD_LIBRARY_PATH
fi

unset SETUPDIR
