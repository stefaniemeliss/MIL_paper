#!/bin/bash
source ~/.bashrc

################################################################################
# extract ROI data for ROI-to-ROI ISFC & create maps for seed-based ISFC
################################################################################ 

# define path
path="/storage/shared/research/cinn/2018/MAGMOT"


# change directory to BIDS folder
BIDS_dir="$path"/rawdata
cd $BIDS_dir

# declare deriv_dir
deriv_dir="$path"/derivatives

# declare ISC dir
ISFC_root="$deriv_dir"/analysis/task_paper/IFSC
mkdir -p $ISFC_root
HPC_dir="$ISFC_root"/aHPC
mkdir -p $HPC_dir
VTA_dir="$ISFC_root"/VTA
mkdir -p $VTA_dir

# declare concat dir
concat_dir="$deriv_dir"/concat

# declare ROI dir
ROI_dir="$deriv_dir"/ROI_masks/output

# define ROI masks
HPC=$ROI_dir/MNI_res-epi_label-aHPC_mask.nii.gz
VTA=$ROI_dir/MNI_res-epi_label-VTASN_mask.nii.gz


# define search and replace string
searchstring="desc-finconcat_bold.nii.gz"
replacestring="maskave_"

# define subjects based on folder names in the BIDS directory
subjects=($(ls -d sub*))
#subjects=(sub-control001 sub-control002 sub-control003) # script development

# get length of subjects
num_subjects=${#subjects[@]}

# define task
task=magictrickwatching

for (( s=0; s<${num_subjects}; s++));
    do
    # define variable subj_id, change directory to where the bold data of subj_id is saved and define file for subj_id
    subj_id=${subjects[$s]}

    cd $concat_dir
    s1_file=($(ls "$subj_id"_task-"$task"_desc-finconcat_bold.nii.gz))

    # define prefix for mask average
	maskave_prefix="${s1_file/$searchstring/$replacestring}"
	VTA_s1="$maskave_prefix"VTA.txt
	HPC_s1="$maskave_prefix"aHPC.txt

    # extract ROI data for subject 1
    if [[ "$subj_id" == *"sub-control001"* ]]; then

		    # extract ROI average time course: for each volume, the average of voxel is computed in saved in txt file
		    3dmaskave -quiet -mask $VTA $s1_file > $ISFC_root/$VTA_s1 # time course of VTA
		    3dmaskave -quiet -mask $HPC $s1_file > $ISFC_root/$HPC_s1 # time course of anterior HPC

    fi

    for (( t=s+1; t<${num_subjects}; t++));
        do

        # define variable subj_corr, change directory to where the bold data of subj_corr is saved and define file for subj_corr
        subj_corr=${subjects[$t]}

        cd $concat_dir
        s2_file=($(ls "$subj_corr"_task-"$task"_desc-finconcat_bold.nii.gz))


        # define prefix
        maskave_prefix="${s2_file/$searchstring/$replacestring}"
        VTA_s2="$maskave_prefix"VTA.txt
        HPC_s2="$maskave_prefix"aHPC.txt

        # extract ROI data for all other subjects when processing subject 1
        if [[ "$subj_id" == *"sub-control001"* ]]; then

		        # extract ROI average time course: for each volume, the average of voxel is computed in saved in txt file
		        3dmaskave -quiet -mask $VTA $s2_file > $ISFC_root/$VTA_s2 # time course of VTA
		        3dmaskave -quiet -mask $HPC $s2_file > $ISFC_root/$HPC_s2 # time course of anterior HPC

        fi

        # determine the prefix for all files
        ISFC_s1_HPC="ISFC_""$subj_id""$subj_corr""_task-"$task"_mask-"$subj_id"_aHPC_z.nii.gz"
        ISFC_s1_VTA="ISFC_""$subj_id""$subj_corr""_task-"$task"_mask-"$subj_id"_VTA_z.nii.gz"
        ISFC_s2_HPC="ISFC_""$subj_id""$subj_corr""_task-"$task"_mask-"$subj_corr"_aHPC_z.nii.gz"
        ISFC_s2_VTA="ISFC_""$subj_id""$subj_corr""_task-"$task"_mask-"$subj_corr"_VTA_z.nii.gz"

        ISFC_HPC="ISFC_""$subj_id""$subj_corr""_task-"$task"_mask-aHPC_z.nii.gz"
        ISFC_VTA="ISFC_""$subj_id""$subj_corr""_task-"$task"_mask-VTA_z.nii.gz"

        # compute ISFC map
        cd $ISFC_root
        if [ ! -f "$ISFC_s1_aHPC" ]; then
        	echo s1_file $s1_file
        	echo s2_file $s2_file
            echo ""

            # compute seed-based map for each ROI for each subject of the pair #

			# calculate ISFC using pearson, save it in .nii.gz, do Fisher-z transformation
			# Usage: 3dTcorr1D [options] xset y1D
            echo ISFC_prefix $ISFC_s1_HPC
            3dTcorr1D -pearson -Fisher -prefix $ISFC_s1_HPC $concat_dir/$s2_file $HPC_s1

			# calculate ISFC using pearson, save it in .nii.gz, do Fisher-z transformation
			# Usage: 3dTcorr1D [options] xset y1D
            echo ISFC_prefix $ISFC_s1_VTA
            3dTcorr1D -pearson -Fisher -prefix $ISFC_s1_VTA $concat_dir/$s2_file $VTA_s1

			# calculate ISFC using pearson, save it in .nii.gz, do Fisher-z transformation
			# Usage: 3dTcorr1D [options] xset y1D
            echo ISFC_prefix $ISFC_s2_HPC
            3dTcorr1D -pearson -Fisher -prefix $ISFC_s2_HPC $concat_dir/$s1_file $HPC_s2

			# calculate ISFC using pearson, save it in .nii.gz, do Fisher-z transformation
			# Usage: 3dTcorr1D [options] xset y1D
            echo ISFC_prefix $ISFC_s2_VTA
            3dTcorr1D -pearson -Fisher -prefix $ISFC_s2_VTA $concat_dir/$s1_file $VTA_s2

            echo ""

            # compute averae of each seed for each pair #
            echo ISFC_prefix $ISFC_HPC
            3dcalc -a $ISFC_s1_HPC -b $ISFC_s2_HPC -expr '(a+b)/2' -prefix $ISFC_HPC -session $HPC_dir

            echo ISFC_prefix $ISFC_VTA
            3dcalc -a $ISFC_s1_VTA -b $ISFC_s2_VTA -expr '(a+b)/2' -prefix $ISFC_VTA -session $VTA_dir

            echo ""

		fi

	done
done
