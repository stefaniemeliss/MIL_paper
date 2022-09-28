#!/bin/bash

source ~/.bashrc

# this script runs the 3dISC command for the ISC-RSA analyses

# define task
task=magictrickwatching

# define path
path="/storage/shared/research/cinn/2018/MAGMOT"
deriv_dir="$path"/derivatives
ana_root="$deriv_dir"/analysis

CS_path=$ana_root/3dClustSim

ana_path="$ana_root"/task_paper/ISFC

###### define masks to loop through ######

# define ROI masks
mask_dir=$deriv_dir/ROI_masks/output

masks=(GM)

mask_files=(sample_label-gm_mask.nii.gz)

num_masks=${#masks[@]}

# define seed ROIs to loop through
seeds=(aHPC VTA)


for seed in "${seeds[@]}"; do

    # change directory
    cd $ana_path/$seed

    # copy 3dClustSim output to directory
    cp $CS_path/ClustSim_magictrickwatching* ./


    # loop over masks to create a 3dISC for each mask
    for (( m=0; m<${num_masks}; m++)); do

	    # define variables
	    mask=${masks[$m]}

	    # run the model for magictrickwatching (no covariates)
        source "$seed"_3dISC_"$task"_"$mask".sh > out_"$seed"_3dISC_"$task"_"$mask".txt

        # Copy the string 'x' (file) into the dataset(s) giving it the name n
        3drefit -atrstring AFNI_CLUSTSIM_NN1_1sided file:ClustSim_magictrickwatching.NN1_1sided.niml \
	    -atrstring AFNI_CLUSTSIM_NN2_1sided file:ClustSim_magictrickwatching.NN2_1sided.niml \
	    -atrstring AFNI_CLUSTSIM_NN3_1sided file:ClustSim_magictrickwatching.NN3_1sided.niml \
	    -atrstring AFNI_CLUSTSIM_NN1_2sided file:ClustSim_magictrickwatching.NN1_2sided.niml \
	    -atrstring AFNI_CLUSTSIM_NN2_2sided file:ClustSim_magictrickwatching.NN2_2sided.niml \
	    -atrstring AFNI_CLUSTSIM_NN3_2sided file:ClustSim_magictrickwatching.NN3_2sided.niml \
	    -atrstring AFNI_CLUSTSIM_NN1_bisided file:ClustSim_magictrickwatching.NN1_bisided.niml \
	    -atrstring AFNI_CLUSTSIM_NN2_bisided file:ClustSim_magictrickwatching.NN2_bisided.niml \
	    -atrstring AFNI_CLUSTSIM_NN3_bisided file:ClustSim_magictrickwatching.NN3_bisided.niml \
	    ISFC_"$seed"_"$task"_"$mask"+tlrc

        covariates=(uniqueCur_highConf uniqueMem_highConf curBetaFull_highConf curBetaRed_highConf)
        
        covariates=(curBetaFull_highConf curBetaRed_highConf)

        num_covariates=${#covariates[@]}


        # loop over covariates to run 3dISC for each of them
	    for (( v=0; v<${num_covariates}; v++)); do
	        
            # define variables
            cov=${covariates[$v]}
		    int=${interaction[$v]}

	        # run the RSA model for each covariate
	        source "$seed"_3dISC_"$task"_"$cov"_"$mask".sh > out_"$seed"_3dISC_"$task"_"$cov"_"$mask".txt

            # Copy the string 'x' (file) into the dataset(s) giving it the name n
            3drefit -atrstring 'AFNI_CLUSTSIM_NN1_1sided' file:ClustSim_magictrickwatching.NN1_1sided.niml \
	        -atrstring 'AFNI_CLUSTSIM_NN2_1sided' file:ClustSim_magictrickwatching.NN2_1sided.niml \
	        -atrstring 'AFNI_CLUSTSIM_NN3_1sided' file:ClustSim_magictrickwatching.NN3_1sided.niml \
	        -atrstring 'AFNI_CLUSTSIM_NN1_2sided' file:ClustSim_magictrickwatching.NN1_2sided.niml \
	        -atrstring 'AFNI_CLUSTSIM_NN2_2sided' file:ClustSim_magictrickwatching.NN2_2sided.niml \
	        -atrstring 'AFNI_CLUSTSIM_NN3_2sided' file:ClustSim_magictrickwatching.NN3_2sided.niml \
	        -atrstring 'AFNI_CLUSTSIM_NN1_bisided' file:ClustSim_magictrickwatching.NN1_bisided.niml \
	        -atrstring 'AFNI_CLUSTSIM_NN2_bisided' file:ClustSim_magictrickwatching.NN2_bisided.niml \
	        -atrstring 'AFNI_CLUSTSIM_NN3_bisided' file:ClustSim_magictrickwatching.NN3_bisided.niml \
	        ISFC_"$seed"_"$task"_"$cov"_"$mask"+tlrc	    

        done
    done

    # rm 3dClustSim output in directory
    rm ./ClustSim_magictrickwatching*
    
done




