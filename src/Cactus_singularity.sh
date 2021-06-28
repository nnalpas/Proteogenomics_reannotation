#!/bin/bash

# Script to execute the MaxQuant singularity

# Time scripts starts
echo "$0"
echo "Start $(date +"%T %d-%m-%Y")."

# Run the singularity runscript section
echo "singularity run ${SCRIPT_FLAGS}"
singularity exec ${SCRIPT_FLAGS}

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."


