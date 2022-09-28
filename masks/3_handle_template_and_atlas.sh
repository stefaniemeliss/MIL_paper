#!/bin/bash
source ~/.bashrc


path="/storage/shared/research/cinn/2018/MAGMOT"

# define directories
deriv_dir="$path"/derivatives
 

# change directory to ROI folder
ROI_root="$deriv_dir"/ROI_masks
ROI_dir="$ROI_root"/output

cd $ROI_dir

# define which template to use and where to find them
template=MNI152_2009_template_SSW.nii.gz
template_path=`@FindAfniDsetPath $template`

# copy MNI template that has been used to normalise data to ROI_dir
cp $template_path/$template ./$template

# define atlas of  interest
glasser=MNI_Glasser_HCP_v1.0.nii.gz
glasser_path=`@FindAfniDsetPath $glasser`
cp $glasser_path/$glasser ./$glasser # copy template to directory

caez=MNI_caez_ml_18+tlrc # Eickhoff-Zilles atlas macro labels (N27-MNI)
caez_path=`@FindAfniDsetPath $caez`
cp $caez_path/$caez* . # copy template to directory

# change atlas list
@AfniEnv -set AFNI_ATLAS_LIST "CA_ML_18_MNI,MNI_Glasser_HCP_v1.0,Brainnetome_1.0,CA_MPM_22_MNI"
@AfniEnv -set AFNI_ATLAS_COLORS MNI_Glasser_HCP_v1.0 # set default atlas

