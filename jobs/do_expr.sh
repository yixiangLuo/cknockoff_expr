#!/bin/bash

EXPERIMENT=$1

# submit experiment jobs
sbatch "${EXPERIMENT}_job.sh"
