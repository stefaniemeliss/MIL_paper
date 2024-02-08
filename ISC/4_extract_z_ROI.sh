#!/bin/bash

# define path
DIR="/mnt/dell_storage/labs/atsuchiyagaito/MIL"

# define derivatves dir
deriv_dir="$DIR"/derivatives

# change directory to BIDS folder
BIDS_dir="$DIR"/rawdata
cd $BIDS_dir

# define subjects based on folder names in the BIDS directory
subjects=($(ls -d sub*))
#subjects=(sub-control001 sub-control002 sub-control003)

# change directory to ISC folder
ISC_dir="$deriv_dir"/analysis/task_paper/ISC
cd $ISC_dir

# Create CSV header
echo "file,aHPC,Caudate,NAcc,VTASN" > ISC_ROI_results.csv

for subj in "${subjects[@]}"; do
    # Nested loop for matched files
    for file in ISC_${subj}*_z.nii.gz; do
        if [[ -f "$file" ]]; then
            aHPC=$(3dmaskave -mask "$deriv_dir"/ROI_masks/output/MNI_res-epi_label-aHPC_mask.nii.gz $file | awk '{print $1}')
            Caudate=$(3dmaskave -mask "$deriv_dir"/ROI_masks/output/MNI_res-epi_label-Caudate_mask.nii.gz $file | awk '{print $1}')
            NAcc=$(3dmaskave -mask "$deriv_dir"/ROI_masks/output/MNI_res-epi_label-NAcc_mask.nii.gz $file | awk '{print $1}')
            VTASN=$(3dmaskave -mask "$deriv_dir"/ROI_masks/output/MNI_res-epi_label-VTASN_mask.nii.gz $file | awk '{print $1}')
            
            # Append results to CSV
            echo "$file,$aHPC,$Caudate,$NAcc,$VTASN" >> ISC_ROI_results.csv
        else
            echo "File not found for pattern ISC_${subj}*_z.nii.gz"
        fi
    done
done

