#!/bin/bash

source ~/.bashrc

# this script creates the 3dISC command the ISC-RSA analyses

# define task
task=magictrickwatching

# define path
path="/storage/shared/research/cinn/2018/MAGMOT"
deriv_dir="$path"/derivatives

ana_path="$deriv_dir"/analysis/task_paper/ISC
mkdir -p $ana_path

script_path=$path/code/task_paper/ISC


# define dataTable
data_table=dataTable_ISC_magictrickwatching.txt

# change directory and copy dataTable
cd $ana_path
cp $script_path/$data_table ./$data_table

###### define masks to loop through ######

# define ROI masks
mask_dir=$deriv_dir/ROI_masks/output

masks=(GM aHPC VTA Caudate NAcc)
mask_files=(sample_label-gm_mask.nii.gz MNI_res-epi_label-aHPC_mask.nii.gz MNI_res-epi_label-VTASN_mask.nii.gz MNI_res-epi_label-Caudate_mask.nii.gz MNI_res-epi_label-NAcc_mask.nii.gz)

num_masks=${#masks[@]}


# loop over masks to create a 3dISC for each mask
for (( m=0; m<${num_masks}; m++)); do

	# define variables
	mask=${masks[$m]}
	mask_f=${mask_files[$m]}

	# specify 3dISC command
	printf "3dISC -prefix ISC_"$task"_"$mask" -jobs 2 \\" > ./3dISC_"$task"_"$mask".sh #name and number of jobs
	printf "\n\t -model 'grp+(1|Subj1)+(1|Subj2)' \\" >> ./3dISC_"$task"_"$mask".sh # model
	printf "\n\t -qVars   'grp' \\" >> ./3dISC_"$task"_"$mask".sh # qVars
	printf "\n\t -gltCode ave '1 0' \\" >> ./3dISC_"$task"_"$mask".sh # contrasts
	printf "\n\t -gltCode G11vG22 '0 1' \\" >> ./3dISC_"$task"_"$mask".sh # contrasts
	printf "\n\t -gltCode G11 '1 0.5' \\" >> ./3dISC_"$task"_"$mask".sh # contrasts
	printf "\n\t -gltCode G22 '1 -0.5' \\" >> ./3dISC_"$task"_"$mask".sh # contrasts
	printf "\n\t -mask $mask_dir/$mask_f \\" >> ./3dISC_"$task"_"$mask".sh # mask
	printf "\n\t -dataTable @$data_table" >> ./3dISC_"$task"_"$mask".sh # datatable

    # define list of covariates and interaction terms
    covariates=(uniqueCur_highConf uniqueMem_highConf curBetaFull_highConf curBetaRed_highConf)
    interaction=(grUniqueCur_highConf grUniqueMem_highConf grCurBetaFull_highConf grCurBetaRed_highConf)
    
    covariates=(curBetaFull_highConf curBetaRed_highConf)
    interaction=(grCurBetaFull_highConf grCurBetaRed_highConf)

    num_covariates=${#covariates[@]}

    # loop over covariates to create a 3dISC for each of them
	for (( v=0; v<${num_covariates}; v++)); do
	    
        # define variables
        cov=${covariates[$v]}
        int=${interaction[$v]}

        # specify 3dISC command including group
        printf "3dISC -prefix ISC_"$task"_"$cov"_"$mask" -jobs 2 \\" > ./3dISC_"$task"_"$cov"_"$mask".sh #name and number of jobs
        printf "\n\t -model 'grp+$cov+$int+(1|Subj1)+(1|Subj2)' \\" >> ./3dISC_"$task"_"$cov"_"$mask".sh # model
        printf "\n\t -qVars   'grp,$cov,$int' \\" >> ./3dISC_"$task"_"$cov"_"$mask".sh # qVars
        printf "\n\t -qVarCenters   '0,0,0' \\" >> ./3dISC_"$task"_"$cov"_"$mask".sh # qVarCenters VALUE, VALUE, VALUE --> this effectively turns off any centering
        printf "\n\t -gltCode ave '1 0 0 0' \\" >> ./3dISC_"$task"_"$cov"_"$mask".sh # contrasts
        printf "\n\t -gltCode G11vG22 '0 1 0 0' \\" >> ./3dISC_"$task"_"$cov"_"$mask".sh # contrasts
        printf "\n\t -gltCode G11 '1 0.5 0 0' \\" >> ./3dISC_"$task"_"$cov"_"$mask".sh # contrasts
        printf "\n\t -gltCode G22 '1 -0.5 0 0' \\" >> ./3dISC_"$task"_"$cov"_"$mask".sh # contrasts
        printf "\n\t -gltCode $cov '0 0 1 0'  \\" >> ./3dISC_"$task"_"$cov"_"$mask".sh # contrasts
        printf "\n\t -gltCode "$cov"G11v"$cov"G22 '0 0 0 1'\\" >> ./3dISC_"$task"_"$cov"_"$mask".sh # contrasts
        printf "\n\t -gltCode "$cov"G11 '0 0 1 0.5'  \\" >> ./3dISC_"$task"_"$cov"_"$mask".sh # contrasts
        printf "\n\t -gltCode "$cov"G22 '0 0 1 -0.5'  \\" >> ./3dISC_"$task"_"$cov"_"$mask".sh # contrasts
        printf "\n\t -mask $mask_dir/$mask_f \\" >> ./3dISC_"$task"_"$cov"_"$mask".sh # mask
        printf "\n\t -dataTable @$data_table" >> ./3dISC_"$task"_"$cov"_"$mask".sh # datatable
        
    done

done






