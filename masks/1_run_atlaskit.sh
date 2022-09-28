#!/bin/bash
source ~/.bashrc

# load module
module load anaconda

# define path
path="/storage/shared/research/cinn/2018/MAGMOT"

# define directories
deriv_dir="$path"/derivatives

# change directory to ROI folder
ROI_root="$deriv_dir"/ROI_masks
input_dir="$ROI_root"/input
anat_dir="$ROI_root"/anatomical_masks
mkdir -p $anat_dir
cd $ROI_root

###### REWARD NETWORK ######

# create masks using atlaskit (https:/github.com/jmtyszka/atlaskit)
# input is Pauli et al (2018) atlas for subcortical structures downloaded from https://osf.io/w8zq2/
# -t set to 0.15 (probability)
# number in the end relates to numerical label

# define strings for output
prefix=MNI_label-
suffix=_mask.nii.gz

# labels for subcortical structures
# 0   Pu	# 1   Ca	# 2   NAC	# 3   EXA	# 4   GPe	# 5   GPi	# 6   SNc
# 7   RN	# 8   SNr	# 9   PBP	# 10  VTA	# 11  VeP	# 12  HN	# 13  HTH	# 14  MN	# 15  STH

# NAcc mask
python /storage/shared/research/cinn/2018/MAGMOT/software/atlaskit-master/create_mask.py -i $input_dir/CIT168toMNI152-2009c_prob.nii.gz -o $anat_dir/"$prefix"NAcc"$suffix" -t 0.15 2

# Caudate mask
python /storage/shared/research/cinn/2018/MAGMOT/software/atlaskit-master/create_mask.py -i $input_dir/CIT168toMNI152-2009c_prob.nii.gz -o $anat_dir/"$prefix"Caudate"$suffix" -t 0.15 1

# SNc mask
python /storage/shared/research/cinn/2018/MAGMOT/software/atlaskit-master/create_mask.py -i $input_dir/CIT168toMNI152-2009c_prob.nii.gz -o $anat_dir/"$prefix"SNc"$suffix" -t 0.15 6

# SNr mask
python /storage/shared/research/cinn/2018/MAGMOT/software/atlaskit-master/create_mask.py -i $input_dir/CIT168toMNI152-2009c_prob.nii.gz -o $anat_dir/"$prefix"SNr"$suffix" -t 0.15 8

# VTA mask
python /storage/shared/research/cinn/2018/MAGMOT/software/atlaskit-master/create_mask.py -i $input_dir/CIT168toMNI152-2009c_prob.nii.gz -o $anat_dir/"$prefix"VTA"$suffix" -t 0.15 10


