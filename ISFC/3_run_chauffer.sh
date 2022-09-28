#!/bin/bash

source ~/.bashrc

# change afni environment vars
@AfniEnv -set AFNI_ATLAS_COLORS MNI_Glasser_HCP_v1.0
@AfniEnv -set AFNI_IMAGE_LABEL_MODE 5
@AfniEnv -set AFNI_VALUE_LABEL YES
@AfniEnv -set AFNI_graph_ggap 0
@AfniEnv -set AFNI_ATLAS_LIST "CA_ML_18_MNI,MNI_Glasser_HCP_v1.0,Brainnetome_1.0,CA_MPM_22_MNI"

# define task
task=magictrickwatching

# define path
path="/storage/shared/research/cinn/2018/MAGMOT"
deriv_dir="$path"/derivatives
ana_root="$deriv_dir"/analysis

ana_path="$ana_root"/task_paper/ISFC

# define ROI masks
mask_dir=$deriv_dir/ROI_masks/output
masks=(GM)
mask_files=(sample_label-gm_mask.nii.gz)
num_masks=${#masks[@]}

# define template
template=MNI152_2009_template_SSW.nii.gz
template_path=`@FindAfniDsetPath $template`
cp $template_path/$template ./$template # copy template to directory

# define yeo networks
yeo=AFNI_Yeo2011_7Networks_MNI152_liberal.nii.gz
yeo_path=$deriv_dir/ROI_masks/input/Yeo_JNeurophysiol11_MNI152
cp $yeo_path/$yeo ./$yeo # copy yeo networks to directory


# define seed ROIs to loop through
seeds=(aHPC VTA)

# define covariates
covariates=(uniqueCur_highConf uniqueMem_highConf curBetaFull_highConf curBetaRed_highConf)
covariates_new=(curiosity memory fullCMLE redCMLE)
#covariates=(uniqueCur_highConf)
#covariates=()
num_covariates=${#covariates[@]}

# define p and k threshold
p_val=0.001
k=20

for seed in "${seeds[@]}"; do

    # change directory
    cd $ana_path/$seed

    # loop over masks to create a 3dISFC for each mask
    for (( m=0; m<${num_masks}; m++)); do

        # define variables
	    mask=${masks[$m]}
	    mask_f=${mask_files[$m]}
	    
	    # copy mask
	    cp $mask_dir/$mask_f ./$mask_f	    

	    ############### save magic trick watching ISFC and reward effect therein ###############

        # define prefixes for ouput files
        clust_ave=clust_ISFC_"$seed"_"$task"_mask-"$mask"_effect-ave.nii.gz
        clust_grp=clust_ISFC_"$seed"_"$task"_mask-"$mask"_effect-grp.nii.gz

        mask_ave=mask_ISFC_"$seed"_"$task"_mask-"$mask"_effect-ave.nii.gz
        mask_grp=mask_ISFC_"$seed"_"$task"_mask-"$mask"_effect-grp.nii.gz
        
        tval_ave=report_ISFC_"$seed"_"$task"_mask-"$mask"_effect-ave_type-tval.txt
        tval_grp=report_ISFC_"$seed"_"$task"_mask-"$mask"_effect-grp_type-tval.txt
    
        effsize_ave=report_ISFC_"$seed"_"$task"_mask-"$mask"_effect-ave_type-effsize.txt
        effsize_grp=report_ISFC_"$seed"_"$task"_mask-"$mask"_effect-grp_type-effsize.txt
        effsize_g11=report_ISFC_"$seed"_"$task"_mask-"$mask"_effect-G11_type-effsize.txt    
        effsize_g22=report_ISFC_"$seed"_"$task"_mask-"$mask"_effect-G22_type-effsize.txt     

        cm_ave=report_ISFC_"$seed"_"$task"_mask-"$mask"_effect-ave_type-HCP_CM.txt # CM = Center of Mass
        min_ave=report_ISFC_"$seed"_"$task"_mask-"$mask"_effect-ave_type-HCP_min.txt # min = Bounding box for the cluster min value
        max_ave=report_ISFC_"$seed"_"$task"_mask-"$mask"_effect-ave_type-HCP_max.txt # max = Bounding box for the cluster max value   
        mi_ave=report_ISFC_"$seed"_"$task"_mask-"$mask"_effect-ave_type-HCP_MI.txt # MI = Maximum Intensity
        ml_hcp_ave=report_ISFC_"$seed"_"$task"_mask-"$mask"_effect-ave_type-HCP_ML.txt # ML = Macro Label    
        ml_caez_ave=report_ISFC_"$seed"_"$task"_mask-"$mask"_effect-ave_type-CAEZ_ML.txt # ML = Macro Label
                
        cm_grp=report_ISFC_"$seed"_"$task"_mask-"$mask"_effect-grp_type-HCP_CM.txt # CM = Center of Mass
        min_grp=report_ISFC_"$seed"_"$task"_mask-"$mask"_effect-grp_type-HCP_min.txt # min = Bounding box for the cluster min value
        max_grp=report_ISFC_"$seed"_"$task"_mask-"$mask"_effect-grp_type-HCP_max.txt # max = Bounding box for the cluster max value   
        mi_grp=report_ISFC_"$seed"_"$task"_mask-"$mask"_effect-grp_type-HCP_MI.txt # MI = Maximum Intensity
        ml_hcp_grp=report_ISFC_"$seed"_"$task"_mask-"$mask"_effect-grp_type-HCP_ML.txt # ML = Macro Label    
        ml_caez_grp=report_ISFC_"$seed"_"$task"_mask-"$mask"_effect-grp_type-CAEZ_ML.txt # ML = Macro Label
        
        # define in and out file prefix
        in=ISFC_"$seed"_"$task"_"$mask"
        out=out_ISFC_"$seed"_"$task"   

	    # threshold output so that only voxel with significant correlation after bonferroni correction survive  --> creates mask
	    3dClusterize -inset "$in"+tlrc -ithr 1 -NN 1 -binary -bisided p=$p_val -clust_nvox $k -pref_map $mask_ave -mask $mask_f > $tval_ave # magictrickwatching
	    3dClusterize -inset "$in"+tlrc -ithr 3 -NN 1 -binary -bisided p=$p_val -clust_nvox $k -pref_map $mask_grp -mask $mask_f > $tval_grp # group effects in magictrick watching
	    # threshold output so that only voxel with significant correlation after bonferroni correction survive  --> creates nary mask for report overlap
	    3dClusterize -inset "$in"+tlrc -idat 0 -ithr 1 -NN 1 -bisided p=$p_val -clust_nvox $k -pref_map $clust_ave -mask $mask_f > $effsize_ave # magictrickwatching
	    3dClusterize -inset "$in"+tlrc -idat 2 -ithr 3 -NN 1 -bisided p=$p_val -clust_nvox $k -pref_map $clust_grp -mask $mask_f > $effsize_grp # group effects in magictrick watching
	    
	    # run report through HCP atlas: ISFC
        if [ -f "$clust_ave" ]; then
	        whereami -tab -coord_file "$effsize_ave"'[1,2,3]' -atlas MNI_Glasser_HCP_v1.0 > $cm_ave # CM = Center of Mass
	        whereami -tab -coord_file "$effsize_ave"'[4,5,6]' -atlas MNI_Glasser_HCP_v1.0 > $min_ave # min = Bounding box for the cluster min value
	        whereami -tab -coord_file "$effsize_ave"'[7,8,9]' -atlas MNI_Glasser_HCP_v1.0 > $max_ave # max = Bounding box for the cluster max value
	        whereami -tab -coord_file "$effsize_ave"'[13,14,15]' -atlas MNI_Glasser_HCP_v1.0 > $mi_ave # MI = Maximum Intensity
	        
	        whereami -omask $clust_ave -atlas MNI_Glasser_HCP_v1.0 > $ml_hcp_ave # ML = macro label
	        whereami -omask $clust_ave -atlas CA_ML_18_MNI > $ml_caez_ave # ML = macro label
	    else
	    # if no clusters are found, create a file only containing zeros
            3dcalc -a $mask_f -expr 'amongst(a,2)' -prefix $mask_ave -datum short
        fi
	
	    # run report through HCP atlas: reward effect
        if [ -f $clust_grp ]; then
            whereami -tab -coord_file "$effsize_grp"'[1,2,3]' -atlas MNI_Glasser_HCP_v1.0 > $cm_grp # CM = Center of Mass
	        whereami -tab -coord_file "$effsize_grp"'[4,5,6]' -atlas MNI_Glasser_HCP_v1.0 > $min_grp # min = Bounding box for the cluster min value
	        whereami -tab -coord_file "$effsize_grp"'[7,8,9]' -atlas MNI_Glasser_HCP_v1.0 > $max_grp # max = Bounding box for the cluster max value
	        whereami -tab -coord_file "$effsize_grp"'[13,14,15]' -atlas MNI_Glasser_HCP_v1.0 > $mi_grp # MI = Maximum Intensity
	        
	        whereami -omask $clust_grp -atlas MNI_Glasser_HCP_v1.0 > $ml_hcp_grp # ML = macro label
	        whereami -omask $clust_grp -atlas CA_ML_18_MNI > $ml_caez_grp # ML = macro label
	        
	        # extract values for each group: use same data set to threshold, but use different data values -idat 
	        # control G11
	        3dClusterize -inset "$in"+tlrc -ithr 3 -idat 4 -NN 1 -bisided p=$p_val -clust_nvox $k -mask $mask_f > $effsize_g11 # group effects in ISFC
	        # experimental G22
	        3dClusterize -inset "$in"+tlrc -ithr 3 -idat 6 -NN 1 -bisided p=$p_val -clust_nvox $k -mask $mask_f > $effsize_g22 # group effects in ISFC

	    else
	    # if no clusters are found, create a file only containing zeros
            3dcalc -a $mask_f -expr 'amongst(a,2)' -prefix $mask_grp -datum short
        fi

	    # use mask to multiply it with output of 3dISFC
	    3dcalc -a $mask_ave -b "$in"+tlrc[0] -expr 'a*b' -prefix masked_ISFC_"$seed".nii.gz -datum float # magictrickwatching
	    3dcalc -a $mask_ave -b "$in"+tlrc[1] -expr 'a*b' -prefix masked_ISFC_"$seed"_t.nii.gz -datum float # magictrickwatching
	    3dcalc -a $mask_grp -b "$in"+tlrc[2] -expr 'a*b' -prefix masked_ISFC_"$seed"_reward.nii.gz -datum float # group effects in magictrickwatching
	    3dcalc -a $mask_grp -b "$in"+tlrc[3] -expr 'a*b' -prefix masked_ISFC_"$seed"_reward_t.nii.gz -datum float # group effects in magictrickwatching

	    # put them all together
	    3dcopy $mask_grp "$out" # group effects in magictrickwatching
	    3dbucket masked_ISFC_"$seed"_reward_t.nii.gz -glueto "$out"+tlrc # group effects in magictrickwatching
	    3dbucket masked_ISFC_"$seed"_reward.nii.gz -glueto "$out"+tlrc # group effects in magictrickwatching
		3dbucket $mask_ave -glueto "$out"+tlrc # covariate mask
	    3dbucket masked_ISFC_"$seed"_t.nii.gz -glueto "$out"+tlrc # magictrickwatching
	    3dbucket masked_ISFC_"$seed".nii.gz -glueto "$out"+tlrc # magictrickwatching

	    # declare sub-brick to be t values
	    3drefit -substatpar 1 fitt  "$out"+tlrc
	    3drefit -substatpar 4 fitt  "$out"+tlrc

	    # rename sub-bricks
	    3drefit -sublabel 0 "ISFC_"$seed -sublabel 1 "ISFC_"$seed"_t" -sublabel 2 "mask_ISFC_"$seed -sublabel 3 "reward" -sublabel 4 "reward_t" -sublabel 5 "mask_reward" "$out"+tlrc

        # results cluster extend thresholded #
        
	    # ISFC map thresholded and cluster extent corrected 
	    fig=fig_ISFC_"$seed"_"$task"_effect-ave
	    func_range=0.05
	    op=7
	    @chauffeur_afni -ulay $template -olay "$out"+tlrc -thr_olay_p2stat $p_val -thr_olay_pside 2sided -set_subbricks 0 0 1 -cbar Reds_and_Blues_Inv -prefix $fig -pbar_saveim "$fig"_pbar -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY -func_range $func_range 

	    # ISFC map thresholded and cluster extent corrected - REWARD EFFECT
	    fig=fig_ISFC_"$seed"_"$task"_effect-grp
	    func_range=0.025
	    op=8
	    @chauffeur_afni -ulay $template -olay "$out"+tlrc -thr_olay_p2stat $p_val -thr_olay_pside 2sided -set_subbricks 0 3 4 -cbar Reds_and_Blues_Inv -prefix $fig -pbar_saveim "$fig"_pbar -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY -func_range $func_range 

        # results with lenient threshold # not used

	    # ISFC map thresholded and cluster extent corrected 
	    #@chauffeur_afni -ulay $template -olay "$in"+tlrc -thr_olay_p2stat 0.2 -thr_olay_pside 2sided -olay_alpha Yes -olay_boxed Yes -set_subbricks 0 0 1 -func_range 0.05 -cbar Reds_and_Blues_Inv -prefix lenient_ISFC_"$seed" -pbar_saveim lenient_ISFC_"$seed"_pbar.png -pbar_dim 64x1351H -opacity 7 -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY

	    # ISFC map thresholded and cluster extent corrected - REWARD EFFECT
	    #@chauffeur_afni -ulay $template -olay "$in"+tlrc -thr_olay_p2stat 0.05 -thr_olay_pside 2sided -olay_alpha Yes -olay_boxed Yes -set_subbricks 0 2 3 -func_range 0.05 -cbar Reds_and_Blues_Inv -prefix lenient_ISFC_"$seed"_reward -pbar_saveim lenient_ISFC_"$seed"_reward_pbar.png -pbar_dim 64x1351H -opacity 7 -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY
	    
	    ################# save main and interaction effects for the covariates in ISRSA analyses #######################


        # loop over covariates to run 3dISFC for each of them
	    for (( v=0; v<${num_covariates}; v++)); do

            # define variables
            cov=${covariates[$v]}
            cov_n=${covariates_new[$v]}
            
            # define prefixes for ouput files
            clust_ave=clust_ISFC_"$seed"_"$cov_n"_mask-"$mask"_effect-main.nii.gz
            clust_grp=clust_ISFC_"$seed"_"$cov_n"_mask-"$mask"_effect-int.nii.gz

            mask_ave=mask_ISFC_"$seed"_"$cov_n"_mask-"$mask"_effect-main.nii.gz
            mask_grp=mask_ISFC_"$seed"_"$cov_n"_mask-"$mask"_effect-int.nii.gz
        
            tval_ave=report_ISFC_"$seed"_"$cov_n"_mask-"$mask"_effect-ave_type-tval.txt
            tval_grp=report_ISFC_"$seed"_"$cov_n"_mask-"$mask"_effect-grp_type-tval.txt
        
            effsize_ave=report_ISFC_"$seed"_"$cov_n"_mask-"$mask"_effect-ave_type-effsize.txt
            effsize_grp=report_ISFC_"$seed"_"$cov_n"_mask-"$mask"_effect-grp_type-effsize.txt
            effsize_g11=report_ISFC_"$seed"_"$cov_n"_mask-"$mask"_effect-G11_type-effsize.txt    
            effsize_g22=report_ISFC_"$seed"_"$cov_n"_mask-"$mask"_effect-G22_type-effsize.txt     

            cm_ave=report_ISFC_"$seed"_"$cov_n"_mask-"$mask"_effect-main_type-HCP_CM.txt # CM = Center of Mass
            min_ave=report_ISFC_"$seed"_"$cov_n"_mask-"$mask"_effect-main_type-HCP_min.txt # min = Bounding box for the cluster min value
            max_ave=report_ISFC_"$seed"_"$cov_n"_mask-"$mask"_effect-main_type-HCP_max.txt # max = Bounding box for the cluster max value   
            mi_ave=report_ISFC_"$seed"_"$cov_n"_mask-"$mask"_effect-main_type-HCP_MI.txt # MI = Maximum Intensity
            ml_hcp_ave=report_ISFC_"$seed"_"$cov_n"_mask-"$mask"_effect-main_type-HCP_ML.txt # ML = Macro Label    
            ml_caez_ave=report_ISFC_"$seed"_"$cov_n"_mask-"$mask"_effect-main_type-CAEZ_ML.txt # ML = Macro Label
                
            cm_grp=report_ISFC_"$seed"_"$cov_n"_mask-"$mask"_effect-int_type-HCP_CM.txt # CM = Center of Mass
            min_grp=report_ISFC_"$seed"_"$cov_n"_mask-"$mask"_effect-int_type-HCP_min.txt # min = Bounding box for the cluster min value
            max_grp=report_ISFC_"$seed"_"$cov_n"_mask-"$mask"_effect-int_type-HCP_max.txt # max = Bounding box for the cluster max value   
            mi_grp=report_ISFC_"$seed"_"$cov_n"_mask-"$mask"_effect-int_type-HCP_MI.txt # MI = Maximum Intensity
            ml_hcp_grp=report_ISFC_"$seed"_"$cov_n"_mask-"$mask"_effect-int_type-HCP_ML.txt # ML = Macro Label    
            ml_caez_grp=report_ISFC_"$seed"_"$cov_n"_mask-"$mask"_effect-int_type-CAEZ_ML.txt # ML = Macro Label
            
            # define in file prefix
            in=ISFC_"$seed"_"$task"_"$cov"_"$mask"
            out=out_ISFC_"$seed"_"$cov_n"
                    
            # extract unthresholded maps
		    3dClusterize -inset "$in"+tlrc -ithr 9 -NN 1 -binary -bisided p=$p_val -clust_nvox $k -pref_map $mask_ave -mask $mask_f > $tval_ave # covariate
		    3dClusterize -inset "$in"+tlrc -ithr 9 -idat 8 -NN 1 -bisided p=$p_val -clust_nvox $k -pref_map $clust_ave -mask $mask_f > $effsize_ave # covariate

	        if [ -f $clust_ave ]; then
	            whereami -tab -coord_file "$effsize_ave"'[1,2,3]' -atlas MNI_Glasser_HCP_v1.0 > $cm_ave # CM = Center of Mass
	            whereami -tab -coord_file "$effsize_ave"'[4,5,6]' -atlas MNI_Glasser_HCP_v1.0 > $min_ave # min = Bounding box for the cluster min value
	            whereami -tab -coord_file "$effsize_ave"'[7,8,9]' -atlas MNI_Glasser_HCP_v1.0 > $max_ave # max = Bounding box for the cluster max value
	            whereami -tab -coord_file "$effsize_ave"'[13,14,15]' -atlas MNI_Glasser_HCP_v1.0 > $mi_ave # MI = Maximum Intensity
	            
	            whereami -omask $clust_ave -atlas MNI_Glasser_HCP_v1.0 > $ml_hcp_ave # ML = macro label
	            whereami -omask $clust_ave -atlas CA_ML_18_MNI > $ml_caez_ave # ML = macro label
	        else # if no clusters are found, create a file only containing zeros    
	            3dcalc -a $mask_f -expr 'amongst(a,2)' -prefix $mask_ave -datum short
            fi

		    # create mask of interaction effect
		    3dClusterize -inset "$in"+tlrc -ithr 11 -NN 1 -binary -bisided p=$p_val -clust_nvox $k -pref_map $mask_grp -mask $mask_f > $tval_grp # covariate-reward interaction
		    3dClusterize -inset "$in"+tlrc -ithr 11 -idat 10 -NN 1 -bisided p=$p_val -clust_nvox $k -pref_map $clust_grp -mask $mask_f > $effsize_grp # covariate-reward interaction
		        
            if [ -f $clust_grp ]; then
	            # run report through HCP atlas: interaction effect
                whereami -tab -coord_file "$effsize_grp"'[1,2,3]' -atlas MNI_Glasser_HCP_v1.0 > $cm_grp # CM = Center of Mass
	            whereami -tab -coord_file "$effsize_grp"'[4,5,6]' -atlas MNI_Glasser_HCP_v1.0 > $min_grp # min = Bounding box for the cluster min value
	            whereami -tab -coord_file "$effsize_grp"'[7,8,9]' -atlas MNI_Glasser_HCP_v1.0 > $max_grp # max = Bounding box for the cluster max value
	            whereami -tab -coord_file "$effsize_grp"'[13,14,15]' -atlas MNI_Glasser_HCP_v1.0 > $mi_grp # MI = Maximum Intensity
	            
	            whereami -omask $clust_grp -atlas MNI_Glasser_HCP_v1.0 > $ml_hcp_grp # ML = macro label
	            whereami -omask $clust_grp -atlas CA_ML_18_MNI > $ml_caez_grp # ML = macro label    
	            
    	        # extract values for each group: use same data set to threshold, but use different data values -idat 
	            # control G11
	            3dClusterize -inset "$in"+tlrc -ithr 11 -idat 12 -NN 1 -bisided p=$p_val -clust_nvox $k -mask $mask_f > $effsize_g11 # group effects in covariate
	            # experimental G22
	            3dClusterize -inset "$in"+tlrc -ithr 11 -idat 14 -NN 1 -bisided p=$p_val -clust_nvox $k -mask $mask_f > $effsize_g22 # group effects in covariate
	                    	            
	        else # if no clusters are found, create a file only containing zeros
	            3dcalc -a $mask_f -expr 'amongst(a,2)' -prefix $mask_grp -datum short
            fi

		    # use mask to multiply it with output
		    3dcalc -a $mask_ave -b "$in"+tlrc[8] -expr 'a*b' -prefix masked_ISFC_"$seed"_"$cov_n"_main.nii.gz -datum float # cov beta
		    3dcalc -a $mask_ave -b "$in"+tlrc[9] -expr 'a*b' -prefix masked_ISFC_"$seed"_"$cov_n"_main_t.nii.gz -datum float # cov t val
		    3dcalc -a $mask_grp -b "$in"+tlrc[10] -expr 'a*b' -prefix masked_ISFC_"$seed"_"$cov_n"_int.nii.gz -datum float # interaction beta
		    3dcalc -a $mask_grp -b "$in"+tlrc[11] -expr 'a*b' -prefix masked_ISFC_"$seed"_"$cov_n"_int_t.nii.gz -datum float # interaction t val
		    
		    # put them all together
		    3dcopy $mask_grp  "$out" # interaction mask
		    3dbucket masked_ISFC_"$seed"_"$cov_n"_int_t.nii.gz -glueto "$out"+tlrc # interaction t value thresholded
		    3dbucket masked_ISFC_"$seed"_"$cov_n"_int.nii.gz -glueto "$out"+tlrc # interaction beta thresholded
		    3dbucket $mask_ave -glueto "$out"+tlrc # covariate mask
		    3dbucket masked_ISFC_"$seed"_"$cov_n"_main_t.nii.gz -glueto "$out"+tlrc # covariate t value thresholded
		    3dbucket masked_ISFC_"$seed"_"$cov_n"_main.nii.gz -glueto "$out"+tlrc # covariate beta thresholded
		    		
		    # declare sub-brick to be t values
		    3drefit -substatpar 1 fitt  "$out"+tlrc
		    3drefit -substatpar 4 fitt  "$out"+tlrc

		    # rename sub-bricks
	        3drefit -sublabel 0 "ISFC_"$seed"_main_""$cov_n" -sublabel 1 "ISFC_"$seed"_main_t_""$cov_n" -sublabel 2 "mask_ISFC_"$seed"_main_"$cov_n -sublabel 3 "ISFC_"$seed"_int_""$cov_n" -sublabel 4 "ISFC_"$seed"_int_t_""$cov_n" -sublabel 5 "mask_ISFC_"$seed"_int_"$cov_n "$out"+tlrc
		    
		    if [[ "$cov" == *"curBeta"* ]]; then
		        func_range=0.5
		    else
		        func_range=0.05
		    fi
		      
            # results cluster extend thresholded #
            
            # plot main effect
		    fig=fig_ISFC_"$seed"_"$cov_n"_effect-main
	        op=8
	        @chauffeur_afni -ulay $template -olay "$out"+tlrc -thr_olay_p2stat $p_val -thr_olay_pside 2sided -set_subbricks 0 0 1 -cbar Reds_and_Blues_Inv -prefix $fig -pbar_saveim "$fig"_pbar -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY -func_range $func_range 
            # plot interaction effect
		    fig=fig_ISFC_"$seed"_"$cov_n"_effect-int
	        op=8
	        @chauffeur_afni -ulay $template -olay "$out"+tlrc -thr_olay_p2stat $p_val -thr_olay_pside 2sided -set_subbricks 0 3 4 -cbar Reds_and_Blues_Inv -prefix $fig -pbar_saveim "$fig"_pbar -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY -func_range $func_range        

        # results with lenient threshold # not used 
            
            # plot main effect
		    #@chauffeur_afni -ulay $template -olay "$in"+tlrc -thr_olay_p2stat 0.05 -thr_olay_pside 2sided -olay_alpha Yes -olay_boxed Yes -set_subbricks 0 8 9 -cbar Reds_and_Blues_Inv -prefix lenient_ISFC_"$seed"_"$cov"_main -pbar_saveim lenient_ISFC_"$seed"_"$cov"_main_pbar.png -pbar_dim 64x1351H -opacity 9 -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY -func_range $func_range
		    # plot interaction effect
            #@chauffeur_afni -ulay $template -olay "$in"+tlrc -thr_olay_p2stat 0.05 -thr_olay_pside 2sided -olay_alpha Yes -olay_boxed Yes -set_subbricks 0 10 11 -cbar Reds_and_Blues_Inv -prefix lenient_ISFC_"$seed"_"$cov"_int -pbar_saveim lenient_ISFC_"$seed"_"$cov"_int_pbar.png -pbar_dim 64x1351H -opacity 9 -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY -func_range $func_range

	       done

        # remove all ISFC maps and clust stim output from directory
        rm mask* clust*         
	       
    done
    
    #### overlay main effects with Yeo networks

    3dcalc -a $yeo -b out_ISFC_"$seed"_"$task"+tlrc.[0] -expr 'ispositive(b)*a' -prefix yeo_"$seed"_"$task"_pos.nii.gz
    3dcalc -a $yeo -b out_ISFC_"$seed"_"$task"+tlrc.[0] -expr 'isnegative(b)*a' -prefix yeo_"$seed"_"$task"_neg.nii.gz

    3dcalc -a $yeo -b out_ISFC_"$seed"_curiosity+tlrc.[2] -expr 'a*b' -prefix yeo_"$seed"_curiosity.nii.gz
    3dcalc -a $yeo -b out_ISFC_"$seed"_memory+tlrc.[2] -expr 'a*b' -prefix yeo_"$seed"_memory.nii.gz

    3dcalc -a $yeo -b out_ISFC_"$seed"_fullCMLE+tlrc.[0] -expr 'ispositive(b)*a' -prefix yeo_"$seed"_clme_pos.nii.gz
    3dcalc -a $yeo -b out_ISFC_"$seed"_fullCMLE+tlrc.[0] -expr 'isnegative(b)*a' -prefix yeo_"$seed"_clme_neg.nii.gz

    op=8
    func_range=7
    
    fig=fig_yeo_"$seed"_"$task"_neg
    @chauffeur_afni -ulay $template -olay yeo_"$seed"_"$task"_neg.nii.gz -set_subbricks 0 0 0 -func_range $func_range -cbar Color_circle_AJJ -prefix $fig -pbar_saveim "$fig"_pbar -pbar_posonly -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY       
    fig=fig_yeo_"$seed"_"$task"_pos
    @chauffeur_afni -ulay $template -olay yeo_"$seed"_"$task"_pos.nii.gz -set_subbricks 0 0 0 -func_range $func_range -cbar Color_circle_AJJ -prefix $fig -pbar_saveim "$fig"_pbar -pbar_posonly -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY       
    
    
    fig=fig_yeo_"$seed"_cmle_neg
    @chauffeur_afni -ulay $template -olay yeo_"$seed"_clme_neg.nii.gz -set_subbricks 0 0 0 -func_range $func_range -cbar Color_circle_AJJ -prefix $fig -pbar_saveim "$fig"_pbar -pbar_posonly -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY       
    fig=fig_yeo_"$seed"_cmle_pos
    @chauffeur_afni -ulay $template -olay yeo_"$seed"_clme_pos.nii.gz -set_subbricks 0 0 0 -func_range $func_range -cbar Color_circle_AJJ -prefix $fig -pbar_saveim "$fig"_pbar -pbar_posonly -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY       

    fig=fig_yeo_"$seed"_curiosity
    @chauffeur_afni -ulay $template -olay yeo_"$seed"_curiosity.nii.gz -set_subbricks 0 0 0 -func_range $func_range -cbar Color_circle_AJJ -prefix $fig -pbar_saveim "$fig"_pbar -pbar_posonly -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY       

    fig=fig_yeo_"$seed"_memory
    @chauffeur_afni -ulay $template -olay yeo_"$seed"_memory.nii.gz -set_subbricks 0 0 0 -func_range $func_range -cbar Color_circle_AJJ -prefix $fig -pbar_saveim "$fig"_pbar -pbar_posonly -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY       

    
done


