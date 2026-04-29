#!/bin/bash
set -euo pipefail

echo "Running on: $(hostname)"
echo "PWD: $(pwd)"
echo "Args: $@"

SCALE="${1}"
USECALI="${2}"
SMEAR="${3}"
NAMEEXT="${4}"

tar -xzf configres.tgz

chmod +x config/makeKUCMSTimeResolution.obj

./config/makeKUCMSTimeResolution.obj "${SCALE}" "${USECALI}" "${SMEAR}" "${NAMEEXT}"

OUTFILE="res2dPlotsTFileCondor.root"
RENAMED="res2dPlots${NAMEEXT}.root"

if [[ -f "${OUTFILE}" ]]; then
    mv "${OUTFILE}" "${RENAMED}"
else
    echo "ERROR: expected output file ${OUTFILE} not found"
    ls -lh
    exit 1
fi

echo "Finished. Output:"
ls -lh "${RENAMED}"
