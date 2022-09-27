#!/bin/bash
source ~/.bashrc

################################################################################
# extract data from task BOLD series after pre-processing
################################################################################


# define path
DIR="/storage/shared/research/cinn/2018/MAGMOT"

# define derivatves dir
deriv_dir="$DIR"/derivatives

# change directory to BIDS folder
BIDS_dir="$DIR"/rawdata
# change directory to the raw NIFI files
cd $BIDS_dir

# define subjects based on folder names in the BIDS directory
subjects=($(ls -d sub*))
#subjects=(sub-control001)

# define task
task=magictrickwatching

# change directory to conact folder
concat_path=$deriv_dir/concat
mkdir -p $concat_path
cd $concat_path


############### create ROI for the secondary visual cortex ###############

# define which template to use and where to find them
template=MNI152_2009_template_SSW.nii.gz
template_path=`@FindAfniDsetPath $template`

# define prefix
roi=V2
prefix=MNI_label-
suffix=_mask.nii.gz
prefix_epi=MNI_res-epi_label-


# extract left and right V2 from Glasser atlas
whereami -mask_atlas_region MNI_Glasser_HCP_v1.0:L_Second_Visual_Area -prefix "$prefix""$roi"_L"$suffix"
whereami -mask_atlas_region MNI_Glasser_HCP_v1.0:R_Second_Visual_Area -prefix "$prefix""$roi"_R"$suffix"

# combine into one ROI
3dcalc -a "$prefix""$roi"_L"$suffix" -b "$prefix""$roi"_R"$suffix" -expr 'a+b' -prefix "$prefix""$roi""$suffix"

# resample anatomical ROI to EPI MNI grid
3dresample -master sub-control001_task-magictrickwatching_desc-initconcat_bold.nii.gz -input "$prefix""$roi""$suffix" -prefix  "$prefix_epi""$roi""$suffix"

# remove all files not needed
rm *_L* *_R* $prefix*

# get all files in an array
#input=sub-control001/func/*initconcat_bold.nii.gz

# create array containing all init concat files
#array=(`ls sub*init*.nii.gz`)


#array=(`find $deriv_dir/sub-*/func -name "*initconcat_bold.nii.gz"`)
#length_array=${#array[@]}
#echo $length_array

# copy files
#for (( f=0; f<${length_array}; f++));
#    do

	# define file_id
#    file_id=${array[$f]}
#	echo $file_id

	# create basename
#	fb_name=$(basename "$file_id" )
#    echo "$fb_name"

	# copy file
#	cp $deriv_dir/$file_id ./$fb_name

#done


# DEFINE MASK here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mask=/storage/shared/research/cinn/2018/MAGMOT/derivatives_old/ROI_masks_old/V2_resampled_EPI.nii.gz
mask="$prefix_epi""$roi""$suffix"

# for each subject in the subjects array
for subject in "${subjects[@]}"; do

	echo $subject

	subj_file=$subject"_task-magictrickwatching_desc-initconcat_bold.nii.gz"

	############### create mean timecourse file for leave one out procedure ###############

	# rename subject file to temp
	mv $subj_file temp.nii.gz

	# compute mean time course of sample without subject
	mean="mean_task-magictrickwatching_desc-initconcat_"$subject".nii.gz"
	3dMean -prefix $mean sub-*.nii.gz

	# change name of temp back
	mv temp.nii.gz $subj_file

	############### extract ROI data from subject and mean file ###############

	# extract single subject data
	3dmaskdump -noijk -mask $mask $subj_file > V2_subject_$subject.txt

	# extract sample mean data
	3dmaskdump -noijk -mask $mask $mean > V2_sample_mean_$subject.txt

done
