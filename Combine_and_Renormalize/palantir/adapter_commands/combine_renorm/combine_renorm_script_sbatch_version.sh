#!/bin/bash

echo "Job created $(date)"
echo "$@"

echo "session start: $(date)"

source /mnt/nasapps/modules/init/bash 
module purge > /dev/null 
module load singularity/3.7.2 > /dev/null 
export PATH=/mnt/projects/CCBR-Pipelines/bin/:${PATH} 
ts=$(date '+%Y-%m-%d') 
echo "$@"

args=$@

tempdir_location=$8

echo "Temp directory: $tempdir_location"
# Need to modify tmp dir to input arg

singularity exec --bind /mnt/:/mnt,$tempdir_location:/tmp \
                        $3 \
                        bash /mnt/projects/CCBR-Pipelines/pipelines/Combine_Renorm_NIDAP/Singularity_run_script.sh $args;
                        
echo "Combine_renorm step finished at $(date)"


echo "Uploading to NIDAP: $(date)"

echo $5

files_to_be_uploaded_list=("$5/seurat_object.rds" \
                          "$5/combine_renorm_console_log.txt" \
                          "$5/combined_renormalized_result.png");
                          
for files_in_upload in $files_to_be_uploaded_list
do
  if [ -f "$files_in_upload" ]; 
  then
    echo "$files_in_upload exists."
  else
    echo "$files_in_upload does not exist."
    #exit 1
  fi
done
    
names_to_be_uploaded_list=("seurat_object.rds" \
                          "combine_renorm_console_log.txt" \
                          "combined_renormalized_result.png");

NIDAP_key=$6

output_dataset_rid=$7

source /mnt/projects/CCBR-Pipelines/palantir/adapter_commands/utils/Upload_to_NIDAP.sh $NIDAP_key \
                                                                                $output_dataset_rid \
                                                                                $files_to_be_uploaded_list \
                                                                                $names_to_be_uploaded_list;
echo "Combine_renorm step uploaded: $(date)"
