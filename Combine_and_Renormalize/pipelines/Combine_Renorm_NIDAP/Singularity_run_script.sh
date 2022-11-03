#!/bin/bash
input_arg=$@
echo "$@"

source /root/miniconda3/etc/profile.d/conda.sh
echo "Conda loaded"
conda activate single-cell-test-Rbase
echo "Conda activated"

echo "Running FRCE_NIDAP_combine_renorm.R"
Rscript /mnt/projects/CCBR-Pipelines/pipelines/Combine_Renorm_NIDAP/FRCE_NIDAP_combine_renorm.R $input_arg

exit