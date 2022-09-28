#!/bin/bash
source ~/.bashrc


path="/storage/shared/research/cinn/2018/MAGMOT"

# define directories
deriv_dir="$path"/derivatives
afniproc_dir="$deriv_dir"/afniproc
 

# change directory to ROI folder
ROI_root="$deriv_dir"/ROI_masks
ROI_dir="$ROI_root"/output
mkdir $ROI_dir
input_dir="$ROI_root"/input
anat_dir="$ROI_root"/anatomical_masks

cd $ROI_dir

# define which template to use and where to find them
template=MNI152_2009_template_SSW.nii.gz
template_path=`@FindAfniDsetPath $template`

# copy MNI template that has been used to normalise data to ROI_dir
3dcopy $template_path/$template $ROI_dir/$template

############### create average EPI and gray matters masks ###############

#cd "$deriv_dir" # change directory to where pre-processed files are
rm $ROI_dir/sample*

task=rest
task=magictrickwatching
task=both

# define prefix
epi_ave=sample_task-"$task"_label-automask_average.nii.gz
epi_mask=sample_label-dilatedGM_mask.nii.gz

gm_ave=sample_task-"$task"_label-gm_average.nii.gz
gm_mask=sample_label-gm_mask.nii.gz

# define string for the masks that are created during pre-processing
mask_ea=*_task-"$task"_space-MNI152NLin2009cAsym_label-dilatedGM_mask.nii.gz
mask_ea=*_space-MNI152NLin2009cAsym_label-dilatedGM_mask.nii.gz

# create average EPI masks restingstate: Combine individual EPI automasks into a group mask
3dMean -datum float -prefix $ROI_dir/$epi_ave  $deriv_dir/sub*/func/$mask_ea
3dcalc -datum float -prefix $ROI_dir/$epi_mask -a $ROI_dir/$epi_ave -expr 'ispositive(a-0.499)'

# remove average
rm $ROI_dir/$epi_ave

# define string for GM mask in EPI resolution and MNI space (note: this mask is the same for task and rest)
mask_gm=*_task-rest.results/follow_ROI_FSGMe+tlrc.BRIK.gz

# create average EPI masks restingstate: Combine individual EPI automasks into a group mask
3dMean -datum float -prefix $ROI_dir/$gm_ave  $deriv_dir/afni_proc/sub-*/$mask_gm
3dcalc -datum float -prefix $ROI_dir/$gm_mask -a $ROI_dir/$gm_ave -expr 'ispositive(a-0.09)'

# remove average
rm $ROI_dir/$gm_ave

################################ create ROIs ################################

prefix=MNI_label-
suffix=_mask.nii.gz

### REWARD NETWORK ###

# note: this assumes that 1_run_atlaskit was executed already

# combine midbrain masks
3dcalc -a $anat_dir/"$prefix"VTA"$suffix" -b $anat_dir/"$prefix"SNc"$suffix" -c $anat_dir/"$prefix"SNr"$suffix" -expr 'a+b+c' -prefix $anat_dir/"$prefix"VTA_SN"$suffix"

# binarise SN/VTA mask
3dmask_tool -input $anat_dir/"$prefix"VTA_SN"$suffix" -prefix $anat_dir/"$prefix"VTASN"$suffix"

# remove SNc/SNr/VTA
rm $anat_dir/*SNc_*
rm $anat_dir/*SNr_*
rm $anat_dir/*VTA_*

### HIPPOCAMPUS ###

# use AFNI to extract HPC mask from Glasser et al 2016 atlas
whereami -mask_atlas_region MNI_Glasser_HCP_v1.0:L_Hippocampus -prefix $anat_dir/"$prefix"HPC_L"$suffix"
whereami -mask_atlas_region MNI_Glasser_HCP_v1.0:R_Hippocampus -prefix $anat_dir/"$prefix"HPC_R"$suffix"

# combine masks for both hemispheres and delete unilateral masks
3dcalc -a $anat_dir/"$prefix"HPC_L"$suffix" -b $anat_dir/"$prefix"HPC_R"$suffix" -expr 'a+b' -prefix $anat_dir/"$prefix"HPC_both"$suffix"

# reorient data set
3dresample -master $template -input $anat_dir/"$prefix"HPC_both"$suffix" -prefix $anat_dir/"$prefix"HPC"$suffix"

rm $anat_dir/*HPC_L* $anat_dir/*HPC_R* $anat_dir/*HPC_both*

# create anterior and posterior HPC mask based on recommendations of Poppenk et al (2013)

# use MNI y=-21 to determine uncal apex and hence landmask to divide anterior and posterior HPC
# this translates into -A -117 (adding -117 planes of zero at anterior edge due to LPI orientation of file) to create posterior HPC
# and into -A -112 [A-P extent 229 voxels --> 229-117=112] to create anterior HPC

# use zeropad to add zeroes to anterior / posterior part of HPC respectively and to hence create pHPC & aHPC
3dZeropad -A -117 -prefix $anat_dir/"$prefix"pHPC"$suffix" $anat_dir/"$prefix"HPC"$suffix"
3dZeropad -P -112 -prefix $anat_dir/"$prefix"aHPC"$suffix" $anat_dir/"$prefix"HPC"$suffix"


################################# RESAMPLE anatomical ROIs to MNI space and EPI grid #################################

# arrays with ROI names
anat=(NAcc Caudate VTASN HPC aHPC pHPC)
#num_anat=${#anat[@]}

# define prefix
prefix_t1=MNI_res-anat_label-
prefix_epi=MNI_res-epi_label-

# for each mask in the anat array
for roi in "${anat[@]}"; do

	    # resample anatomical ROI to EPI MNI grid
	    3dresample -master "$ROI_dir"/$gm_mask -input $anat_dir/"$prefix""$roi""$suffix" -prefix $ROI_dir/"$prefix_epi""$roi""$suffix"
	
done

3dresample -master "$ROI_dir"/$gm_mask -input $ROI_dir/$template -prefix $ROI_dir/"$prefix_epi"SSW_template.nii.gz





