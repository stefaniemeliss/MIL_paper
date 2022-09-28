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

ana_path="$ana_root"/task_paper/ISC
cd $ana_path    

# define ROI masks
mask_dir=$deriv_dir/ROI_masks/output
masks=(GM aHPC VTA NAcc Caudate)
#masks=(GM)
#masks=(Caudate)
mask_files=(sample_label-gm_mask.nii.gz MNI_res-epi_label-aHPC_mask.nii.gz MNI_res-epi_label-VTASN_mask.nii.gz MNI_res-epi_label-NAcc_mask.nii.gz MNI_res-epi_label-Caudate_mask.nii.gz)
#mask_files=(sample_label-gm_mask.nii.gz)
#mask_files=(MNI_res-epi_label-Caudate_mask.nii.gz)
num_masks=${#masks[@]}

# define template
template=MNI152_2009_template_SSW.nii.gz
template_path=`@FindAfniDsetPath $template`
cp $template_path/$template ./$template # copy template to directory

# define yeo networks
yeo=AFNI_Yeo2011_7Networks_MNI152_liberal.nii.gz
yeo_path=$deriv_dir/ROI_masks/input/Yeo_JNeurophysiol11_MNI152
cp $yeo_path/$yeo ./$yeo # copy yeo networks to directory

# define covariates
covariates=(uniqueCur_highConf uniqueMem_highConf curBetaFull_highConf curBetaRed_highConf)
covariates_new=(curiosity memory fullCMLE redCMLE)
#covariates=(uniqueMem_highConf)
#covariates_new=(memory)
num_covariates=${#covariates[@]}

# loop over masks to create a 3dISC for each mask
for (( m=0; m<${num_masks}; m++)); do

    # define variables
    mask=${masks[$m]}
    mask_f=${mask_files[$m]}

    # copy mask
    cp $mask_dir/$mask_f ./$mask_f
    
    ############### save magic trick watching ISC and reward effect therein ###############
    
    # define prefixes for ouput files
    clust_ave=clust_ISC_"$task"_mask-"$mask"_effect-ave.nii.gz
    clust_grp=clust_ISC_"$task"_mask-"$mask"_effect-grp.nii.gz

    mask_ave=mask_ISC_"$task"_mask-"$mask"_effect-ave.nii.gz
    mask_grp=mask_ISC_"$task"_mask-"$mask"_effect-grp.nii.gz

    tval_ave=report_ISC_"$task"_mask-"$mask"_effect-ave_type-tval.txt
    tval_grp=report_ISC_"$task"_mask-"$mask"_effect-grp_type-tval.txt
    
    effsize_ave=report_ISC_"$task"_mask-"$mask"_effect-ave_type-effsize.txt
    effsize_grp=report_ISC_"$task"_mask-"$mask"_effect-grp_type-effsize.txt
    effsize_g11=report_ISC_"$task"_mask-"$mask"_effect-G11_type-effsize.txt    
    effsize_g22=report_ISC_"$task"_mask-"$mask"_effect-G22_type-effsize.txt     
    
    cm_ave=report_ISC_"$task"_mask-"$mask"_effect-ave_type-HCP_CM.txt # CM = Center of Mass
    min_ave=report_ISC_"$task"_mask-"$mask"_effect-ave_type-HCP_min.txt # min = Bounding box for the cluster min value
    max_ave=report_ISC_"$task"_mask-"$mask"_effect-ave_type-HCP_max.txt # max = Bounding box for the cluster max value   
    mi_ave=report_ISC_"$task"_mask-"$mask"_effect-ave_type-HCP_MI.txt # MI = Maximum Intensity
    ml_hcp_ave=report_ISC_"$task"_mask-"$mask"_effect-ave_type-HCP_ML.txt # ML = Macro Label    
    ml_caez_ave=report_ISC_"$task"_mask-"$mask"_effect-ave_type-CAEZ_ML.txt # ML = Macro Label
            
    cm_grp=report_ISC_"$task"_mask-"$mask"_effect-grp_type-HCP_CM.txt # CM = Center of Mass
    min_grp=report_ISC_"$task"_mask-"$mask"_effect-grp_type-HCP_min.txt # min = Bounding box for the cluster min value
    max_grp=report_ISC_"$task"_mask-"$mask"_effect-grp_type-HCP_max.txt # max = Bounding box for the cluster max value   
    mi_grp=report_ISC_"$task"_mask-"$mask"_effect-grp_type-HCP_MI.txt # MI = Maximum Intensity
    ml_hcp_grp=report_ISC_"$task"_mask-"$mask"_effect-grp_type-HCP_ML.txt # ML = Macro Label    
    ml_caez_grp=report_ISC_"$task"_mask-"$mask"_effect-grp_type-CAEZ_ML.txt # ML = Macro Label
    
    # define in file prefix
    in=ISC_"$task"_"$mask"
    
    # define out file prefix
    out=out_ISC_"$task"_"$mask"

    
    # make adaptions depending on whether this is whole brain or roi based analysis
    if [[ "$mask" == *"GM"* ]]; then
        # define thresholds
        p_val=0.001
        k=20
        # threshold output so that only voxel with significant correlation after bonferroni correction survive  --> creates binary mask to multiply with
	    3dClusterize -inset "$in"+tlrc -ithr 1 -NN 1 -binary -bisided p=$p_val -clust_nvox $k -pref_map $mask_ave -mask $mask_f > $tval_ave # magictrickwatching
	    3dClusterize -inset "$in"+tlrc -ithr 3 -NN 1 -binary -bisided p=$p_val -clust_nvox $k -pref_map $mask_grp -mask $mask_f > $tval_grp # group effects in magictrick watching
	    # threshold output so that only voxel with significant correlation after bonferroni correction survive  --> creates nary mask for report overlap
	    3dClusterize -inset "$in"+tlrc -ithr 1 -idat 0 -NN 1 -bisided p=$p_val -clust_nvox $k -pref_map $clust_ave -mask $mask_f > $effsize_ave # magictrickwatching
	    3dClusterize -inset "$in"+tlrc -ithr 3 -idat 2 -NN 1 -bisided p=$p_val -clust_nvox $k -pref_map $clust_grp -mask $mask_f > $effsize_grp # group effects in magictrick watching
    else
        # define thresholds
        p_val=0.05
        k=5
        # map q values to z values
        3dFDR -input "$in"+tlrc -prefix FDR_"$in"
        # threshold output so that only voxel with FDR(q) < 0.05 survive  --> creates binary mask
        # use the output of the 3dFDR command to threshold data based on z values, but output the effect size measures
	    3dClusterize -inset FDR_"$in"+tlrc -ithr 1 -idat 0 -NN 1 -binary -bisided p=$p_val -clust_nvox $k -pref_map $mask_ave -mask $mask_f > $effsize_ave # magictrickwatching
	    3dClusterize -inset FDR_"$in"+tlrc -ithr 3 -idat 2 -NN 1 -binary -bisided p=$p_val -clust_nvox $k -pref_map $mask_grp -mask $mask_f > $effsize_grp # group effects in magictrick watching

        # use the mask created above to extract t values 
	    3dClusterize -inset "$in"+tlrc -ithr 1 -NN 1 -bisided p=0.99 -mask $mask_ave > $tval_ave # magictrickwatching
	    3dClusterize -inset "$in"+tlrc -ithr 3 -NN 1 -bisided p=0.99 -mask $mask_grp > $tval_grp # group effects in magictrick watching
	    
        # use binary mask to multiply it with ISC data
		3dcalc -a $mask_ave -b "$in"+tlrc.[0] -expr 'a*b' -prefix masked_ISC_"$mask".nii.gz -datum float # magictrickwatching
		3dcalc -a $mask_ave -b "$in"+tlrc.[1] -expr 'a*b' -prefix masked_ISC_t_"$mask".nii.gz -datum float # magictrickwatching
	    
	    # put them all together
	    3dcopy $mask_ave "$out" # mask ISC in ROI
	    3dbucket masked_ISC_t_"$mask".nii.gz -glueto "$out"+tlrc # magictrickwatching
	    3dbucket masked_ISC_"$mask".nii.gz -glueto "$out"+tlrc # magictrickwatching
	    
	    # rename sub-bricks
		3drefit -sublabel 0 "ISC_""$mask" -sublabel 1 "ISC_t_""$mask" -sublabel 2 "mask_ISC""$mask" "$out"+tlrc
		3drefit -substatpar 1 fitt  "$out"+tlrc

    fi
	
	# run report through HCP atlas: magic trick watching
    if [ -f "$clust_ave" ]; then
	    whereami -tab -coord_file "$effsize_ave"'[1,2,3]' -atlas MNI_Glasser_HCP_v1.0 > $cm_ave # CM = Center of Mass
	    whereami -tab -coord_file "$effsize_ave"'[4,5,6]' -atlas MNI_Glasser_HCP_v1.0 > $min_ave # min = Bounding box for the cluster min value
	    whereami -tab -coord_file "$effsize_ave"'[7,8,9]' -atlas MNI_Glasser_HCP_v1.0 > $max_ave # max = Bounding box for the cluster max value
	    whereami -tab -coord_file "$effsize_ave"'[13,14,15]' -atlas MNI_Glasser_HCP_v1.0 > $mi_ave # MI = Maximum Intensity
	    
	    whereami -omask $clust_ave -atlas MNI_Glasser_HCP_v1.0 > $ml_hcp_ave # ML = macro label
	    whereami -omask $clust_ave -atlas CA_ML_18_MNI > $ml_caez_ave # ML = macro label
    fi
	
	# run report through HCP atlas: reward effect
    if [ -f "$clust_grp" ]; then
        whereami -tab -coord_file "$effsize_grp"'[1,2,3]' -atlas MNI_Glasser_HCP_v1.0 > $cm_grp # CM = Center of Mass
	    whereami -tab -coord_file "$effsize_grp"'[4,5,6]' -atlas MNI_Glasser_HCP_v1.0 > $min_grp # min = Bounding box for the cluster min value
	    whereami -tab -coord_file "$effsize_grp"'[7,8,9]' -atlas MNI_Glasser_HCP_v1.0 > $max_grp # max = Bounding box for the cluster max value
	    whereami -tab -coord_file "$effsize_grp"'[13,14,15]' -atlas MNI_Glasser_HCP_v1.0 > $mi_grp # MI = Maximum Intensity
	    
	    whereami -omask $clust_grp -atlas MNI_Glasser_HCP_v1.0 > $ml_hcp_grp # ML = macro label
	    whereami -omask $clust_grp -atlas CA_ML_18_MNI > $ml_caez_grp # ML = macro label
	    
	    # extract values for each group: use same data set to threshold, but use different data values -idat 
	    # control G11
	    3dClusterize -inset "$in"+tlrc -ithr 3 -idat 4 -NN 1 -bisided p=$p_val -clust_nvox $k -mask $mask_f > $effsize_g11 # group effects in magictrick watching
	    # experimental G22
	    3dClusterize -inset "$in"+tlrc -ithr 3 -idat 6 -NN 1 -bisided p=$p_val -clust_nvox $k -mask $mask_f > $effsize_g22 # group effects in magictrick watching
	    
    fi
    
    ############### create plots of whole brain effects ###############    
    
    if [[ "$mask" == *"GM"* ]]; then
    
	    # use mask to multiply it with output of 3dISC
	    3dcalc -a $mask_ave -b "$in"+tlrc.[0] -expr 'a*b' -prefix masked_ISC.nii.gz -datum float # magictrickwatching
	    3dcalc -a $mask_ave -b "$in"+tlrc.[1] -expr 'a*b' -prefix masked_ISC_t.nii.gz -datum float # magictrickwatching
	    3dcalc -a $mask_grp -b "$in"+tlrc.[2] -expr 'a*b' -prefix masked_ISC_reward.nii.gz -datum float # group effects in magictrickwatching
	    3dcalc -a $mask_grp -b "$in"+tlrc.[3] -expr 'a*b' -prefix masked_ISC_reward_t.nii.gz -datum float # group effects in magictrickwatching

	    # put them all together
	    3dcopy $mask_grp "$out" # group effects in magictrickwatching
	    3dbucket masked_ISC_reward_t.nii.gz -glueto "$out"+tlrc # group effects in magictrickwatching
	    3dbucket masked_ISC_reward.nii.gz -glueto "$out"+tlrc # group effects in magictrickwatching
		3dbucket $mask_ave -glueto "$out"+tlrc # covariate mask
	    3dbucket masked_ISC_t.nii.gz -glueto "$out"+tlrc # magictrickwatching
	    3dbucket masked_ISC.nii.gz -glueto "$out"+tlrc # magictrickwatching

	    # declare sub-brick to be t values
	    3drefit -substatpar 1 fitt  "$out"+tlrc
	    3drefit -substatpar 4 fitt  "$out"+tlrc

	    # rename sub-bricks
		3drefit -sublabel 0 "ISC" -sublabel 1 "ISC_t" -sublabel 2 "mask_ISC" -sublabel 3 "reward" -sublabel 4 "reward_t" -sublabel 5 "mask_reward" "$out"+tlrc

        # results cluster extend thresholded #

	    # ISC map thresholded and cluster extent corrected 
	    fig=fig_ISC_"$task"_effect-ave
	    func_range=0.2
	    op=7
	    @chauffeur_afni -ulay $template -olay "$out"+tlrc -thr_olay_p2stat $p_val -thr_olay_pside 2sided -set_subbricks 0 0 1 -cbar Reds_and_Blues_Inv -func_range $func_range -prefix $fig -pbar_saveim "$fig"_pbar -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY

	    # ISC map thresholded and cluster extent corrected - REWARD EFFECT
	    fig=fig_ISC_"$task"_effect-grp
	    func_range=0.05
	    op=8
	    @chauffeur_afni -ulay $template -olay "$out"+tlrc -thr_olay_p2stat $p_val -thr_olay_pside 2sided -set_subbricks 0 3 4 -cbar Reds_and_Blues_Inv -func_range $func_range -prefix $fig -pbar_saveim "$fig"_pbar -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY

        # results with lenient threshold # not used

	    # ISC map thresholded and cluster extent corrected 
	    #@chauffeur_afni -ulay $template -olay ISC_"$task"_GM+tlrc -thr_olay_p2stat 0.2 -thr_olay_pside 2sided -olay_alpha Yes -olay_boxed Yes -set_subbricks 0 0 1 -func_range 0.05 -cbar Reds_and_Blues_Inv -prefix lenient_ISC -pbar_saveim lenient_ISC_pbar.png -pbar_dim 64x1351H -opacity 7 -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY

	    # ISC map thresholded and cluster extent corrected - REWARD EFFECT
	    #@chauffeur_afni -ulay $template -olay ISC_"$task"_GM+tlrc -thr_olay_p2stat 0.05 -thr_olay_pside 2sided -olay_alpha Yes -olay_boxed Yes -set_subbricks 0 2 3 -func_range 0.05 -cbar Reds_and_Blues_Inv -prefix lenient_ISC_reward -pbar_saveim lenient_ISC_reward_pbar.png -pbar_dim 64x1351H -opacity 7 -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY
	
	fi	    
	
	################# save main and interaction effects for the covariates in ISRSA analyses #######################

    # loop over covariates to run 3dISC for each of them
    for (( v=0; v<${num_covariates}; v++)); do

        # define variables
        cov=${covariates[$v]}
        cov_n=${covariates_new[$v]}
        
        # define prefixes for ouput files
        clust_ave=clust_ISC_"$cov_n"_mask-"$mask"_effect-main.nii.gz
        clust_grp=clust_ISC_"$cov_n"_mask-"$mask"_effect-int.nii.gz

        mask_ave=mask_ISC_"$cov_n"_mask-"$mask"_effect-main.nii.gz
        mask_grp=mask_ISC_"$cov_n"_mask-"$mask"_effect-int.nii.gz
    
        tval_ave=report_ISC_"$cov_n"_mask-"$mask"_effect-ave_type-tval.txt
        tval_grp=report_ISC_"$cov_n"_mask-"$mask"_effect-grp_type-tval.txt
    
        effsize_ave=report_ISC_"$cov_n"_mask-"$mask"_effect-ave_type-effsize.txt
        effsize_grp=report_ISC_"$cov_n"_mask-"$mask"_effect-grp_type-effsize.txt
        effsize_g11=report_ISC_"$cov_n"_mask-"$mask"_effect-G11_type-effsize.txt    
        effsize_g22=report_ISC_"$cov_n"_mask-"$mask"_effect-G22_type-effsize.txt     

        cm_ave=report_ISC_"$cov_n"_mask-"$mask"_effect-main_type-HCP_CM.txt # CM = Center of Mass
        min_ave=report_ISC_"$cov_n"_mask-"$mask"_effect-main_type-HCP_min.txt # min = Bounding box for the cluster min value
        max_ave=report_ISC_"$cov_n"_mask-"$mask"_effect-main_type-HCP_max.txt # max = Bounding box for the cluster max value   
        mi_ave=report_ISC_"$cov_n"_mask-"$mask"_effect-main_type-HCP_MI.txt # MI = Maximum Intensity
        ml_hcp_ave=report_ISC_"$cov_n"_mask-"$mask"_effect-main_type-HCP_ML.txt # ML = Macro Label    
        ml_caez_ave=report_ISC_"$cov_n"_mask-"$mask"_effect-main_type-CAEZ_ML.txt # ML = Macro Label
            
        cm_grp=report_ISC_"$cov_n"_mask-"$mask"_effect-int_type-HCP_CM.txt # CM = Center of Mass
        min_grp=report_ISC_"$cov_n"_mask-"$mask"_effect-int_type-HCP_min.txt # min = Bounding box for the cluster min value
        max_grp=report_ISC_"$cov_n"_mask-"$mask"_effect-int_type-HCP_max.txt # max = Bounding box for the cluster max value   
        mi_grp=report_ISC_"$cov_n"_mask-"$mask"_effect-int_type-HCP_MI.txt # MI = Maximum Intensity
        ml_hcp_grp=report_ISC_"$cov_n"_mask-"$mask"_effect-int_type-HCP_ML.txt # ML = Macro Label    
        ml_caez_grp=report_ISC_"$cov_n"_mask-"$mask"_effect-int_type-CAEZ_ML.txt # ML = Macro Label
        
        # define in file prefix
        in=ISC_"$task"_"$cov"_"$mask"
        
        # make adaptions depending on whether this is whole brain or roi based analysis
		if [[ "$mask" == *"GM"* ]]; then
		
		    # threshold output so that only voxel with significant correlation after bonferroni correction survive  --> creates mask
			3dClusterize -inset "$in"+tlrc -ithr 9 -NN 1 -binary -bisided p=$p_val -clust_nvox $k -pref_map $mask_ave -mask $mask_f > $tval_ave # covariate
			3dClusterize -inset "$in"+tlrc -ithr 11 -NN 1 -binary -bisided p=$p_val -clust_nvox $k -pref_map $mask_grp -mask $mask_f > $tval_grp # covariate-reward interaction
			# threshold output so that only voxel with significant correlation after bonferroni correction survive  --> creates nary mask for report overlap
			3dClusterize -inset "$in"+tlrc -ithr 9 -idat 8 -NN 1 -bisided p=$p_val -clust_nvox $k -pref_map $clust_ave -mask $mask_f > $effsize_ave # covariate
			3dClusterize -inset "$in"+tlrc -ithr 11 -idat 10 -NN 1 -bisided p=$p_val -clust_nvox $k -pref_map $clust_grp -mask $mask_f > $effsize_grp # covariate-reward interaction
		else
		
		    # map q values to z values
            3dFDR -input "$in"+tlrc -prefix FDR_"$in"
            # threshold output so that only voxel with FDR(q) < 0.05 survive  --> creates binary mask
            # use the output of the 3dFDR command to threshold data based on z values, but output the effect size measures
	        3dClusterize -inset FDR_"$in"+tlrc -ithr 9 -idat 8 -NN 1 -binary -bisided p=$p_val -clust_nvox $k -pref_map $mask_ave -mask $mask_f > $effsize_ave # magictrickwatching
	        3dClusterize -inset FDR_"$in"+tlrc -ithr 11 -idat 10 -NN 1 -binary -bisided p=$p_val -clust_nvox $k -pref_map $mask_grp -mask $mask_f > $effsize_grp # group effects in magictrick watching

            # use the mask created above to extract t values 
	        3dClusterize -inset "$in"+tlrc -ithr 9 -NN 1 -bisided p=0.99 -mask $mask_ave > $tval_ave # magictrickwatching
	        3dClusterize -inset "$in"+tlrc -ithr 11 -NN 1 -bisided p=0.99 -mask $mask_grp > $tval_grp # group effects in magictrick watching

		fi
		    
	    # run report through HCP atlas: main effect
        if [ -f "$clust_ave" ]; then
	        whereami -tab -coord_file "$effsize_ave"'[1,2,3]' -atlas MNI_Glasser_HCP_v1.0 > $cm_ave # CM = Center of Mass
	        whereami -tab -coord_file "$effsize_ave"'[4,5,6]' -atlas MNI_Glasser_HCP_v1.0 > $min_ave # min = Bounding box for the cluster min value
	        whereami -tab -coord_file "$effsize_ave"'[7,8,9]' -atlas MNI_Glasser_HCP_v1.0 > $max_ave # max = Bounding box for the cluster max value
	        whereami -tab -coord_file "$effsize_ave"'[13,14,15]' -atlas MNI_Glasser_HCP_v1.0 > $mi_ave # MI = Maximum Intensity
	        
	        whereami -omask $clust_ave -atlas MNI_Glasser_HCP_v1.0 > $ml_hcp_ave # ML = macro label
	        whereami -omask $clust_ave -atlas CA_ML_18_MNI > $ml_caez_ave # ML = macro label
        fi
	    
	    # run report through HCP atlas: interaction effect
        if [ -f "$clust_grp" ]; then
            whereami -tab -coord_file "$effsize_grp"'[1,2,3]' -atlas MNI_Glasser_HCP_v1.0 > $cm_grp # CM = Center of Mass
	        whereami -tab -coord_file "$effsize_grp"'[4,5,6]' -atlas MNI_Glasser_HCP_v1.0 > $min_grp # min = Bounding box for the cluster min value
	        whereami -tab -coord_file "$effsize_grp"'[7,8,9]' -atlas MNI_Glasser_HCP_v1.0 > $max_grp # max = Bounding box for the cluster max value
	        whereami -tab -coord_file "$effsize_grp"'[13,14,15]' -atlas MNI_Glasser_HCP_v1.0 > $mi_grp # MI = Maximum Intensity
	        
	        whereami -omask $clust_grp -atlas MNI_Glasser_HCP_v1.0 > $ml_hcp_grp # ML = macro label
	        whereami -omask $clust_grp -atlas CA_ML_18_MNI > $ml_caez_grp # ML = macro label
	        
    	    # extract values for each group: use same data set to threshold, but use different data values -idat 
	        # control G11
	        3dClusterize -inset "$in"+tlrc -ithr 11 -idat 12 -NN 1 -bisided p=$p_val -clust_nvox $k -mask $mask_f > $effsize_g11 # group effects in magictrick watching
	        # experimental G22
	        3dClusterize -inset "$in"+tlrc -ithr 11 -idat 14 -NN 1 -bisided p=$p_val -clust_nvox $k -mask $mask_f > $effsize_g22 # group effects in magictrick watching
	        
        fi
	        
        ############### create plots of whole brain effects ###############    
    
        if [[ "$mask" == *"GM"* ]]; then
        
            # define out file prefix
            out=out_ISC_"$cov_n"_mask-"$mask"

		    # use mask to multiply it with output
		    3dcalc -a $mask_ave -b "$in"+tlrc.[8] -expr 'a*b' -prefix masked_ISC_"$cov_n"_mask-"$mask"_main.nii.gz -datum float # cov beta
		    3dcalc -a $mask_ave -b "$in"+tlrc.[9] -expr 'a*b' -prefix masked_ISC_"$cov_n"_mask-"$mask"_main_t.nii.gz -datum float # cov t val
		    3dcalc -a $mask_grp -b "$in"+tlrc.[10] -expr 'a*b' -prefix masked_ISC_"$cov_n"_mask-"$mask"_int.nii.gz -datum float # interaction beta
		    3dcalc -a $mask_grp -b "$in"+tlrc.[11] -expr 'a*b' -prefix masked_ISC_"$cov_n"_mask-"$mask"_int_t.nii.gz -datum float # interaction t val
		    
		    # put them all together
		    3dcopy $mask_grp  "$out" # interaction mask
		    3dbucket masked_ISC_"$cov_n"_mask-"$mask"_int_t.nii.gz -glueto "$out"+tlrc # interaction t value thresholded
		    3dbucket masked_ISC_"$cov_n"_mask-"$mask"_int.nii.gz -glueto "$out"+tlrc # interaction beta thresholded
		    3dbucket $mask_ave -glueto "$out"+tlrc # covariate mask
		    3dbucket masked_ISC_"$cov_n"_mask-"$mask"_main_t.nii.gz -glueto "$out"+tlrc # covariate t value thresholded
		    3dbucket masked_ISC_"$cov_n"_mask-"$mask"_main.nii.gz -glueto "$out"+tlrc # covariate beta thresholded
				    
		    # declare sub-brick to be t values
		    3drefit -substatpar 1 fitt  "$out"+tlrc
		    3drefit -substatpar 4 fitt  "$out"+tlrc

		    # rename sub-bricks
		    3drefit -sublabel 0 "main_""$cov_n" -sublabel 1 "main_t_""$cov_n" -sublabel 2 "mask_cov" -sublabel 3 "int_""$cov_n" -sublabel 4 "int_t_""$cov_n" -sublabel 5 "mask_int" "$out"+tlrc
		    
		    if [[ "$cov" == *"curBeta"* ]]; then
		        func_range=1
		    else
		        func_range=0.15
		    fi
		      
		    # results cluster extend thresholded # 

		    # plot main effect
		    fig=fig_ISC_"$cov_n"_mask-"$mask"_effect-main
	        op=8
	        @chauffeur_afni -ulay $template -olay "$out"+tlrc -thr_olay_p2stat $p_val -thr_olay_pside 2sided -set_subbricks 0 0 1 -cbar Reds_and_Blues_Inv -func_range $func_range -prefix $fig -pbar_saveim "$fig"_pbar -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY
            # plot interaction effect
		    fig=fig_ISC_"$cov_n"_mask-"$mask"_effect-int
	        op=8
	        @chauffeur_afni -ulay $template -olay "$out"+tlrc -thr_olay_p2stat $p_val -thr_olay_pside 2sided -set_subbricks 0 3 4 -cbar Reds_and_Blues_Inv -func_range $func_range -prefix $fig -pbar_saveim "$fig"_pbar -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY	        

		    # results with lenient threshold # not used
		    
		    # plot main effect
		    #@chauffeur_afni -ulay $template -olay "$in"+tlrc -thr_olay_p2stat 0.05 -thr_olay_pside 2sided -olay_alpha Yes -olay_boxed Yes -set_subbricks 0 8 9 -cbar Reds_and_Blues_Inv -prefix lenient_ISC_"$cov"_main -pbar_saveim lenient_ISC_"$cov"_main_pbar.png -pbar_dim 64x1351H -opacity 9 -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY -func_range $func_range
		    # plot interaction effect
		    #@chauffeur_afni -ulay $template -olay "$in"+tlrc -thr_olay_p2stat 0.05 -thr_olay_pside 2sided -olay_alpha Yes -olay_boxed Yes -set_subbricks 0 10 11 -cbar Reds_and_Blues_Inv -prefix lenient_ISC_"$cov"_int -pbar_saveim lenient_ISC_"$cov"_int_pbar.png -pbar_dim 64x1351H -opacity 9 -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY -func_range $func_range
		    
		fi
		
        ############### create plot for memory and CMLE effects in ROIs ###############    
        
        if [[ "$cov_n" == *"memory"* && "$mask" == *"Caudate"* ]] || [[ "$cov_n" == *"CMLE"* && "$mask" != *"GM"* ]]; then  
            
            # define out file prefix
            out=out_ISC_"$cov_n"_mask-"$mask"

		    # use mask to multiply it with output
		    		    # use mask to multiply it with output
		    3dcalc -a $mask_ave -b "$in"+tlrc.[8] -expr 'a*b' -prefix masked_ISC_"$cov_n"_mask-"$mask"_main.nii.gz -datum float # cov beta
		    3dcalc -a $mask_ave -b "$in"+tlrc.[9] -expr 'a*b' -prefix masked_ISC_"$cov_n"_mask-"$mask"_main_z.nii.gz -datum float # cov z val (after FDR conversion)
		    
		    # put them all together
		    3dcopy $mask_ave  "$out" # interaction mask
		    3dbucket $mask_ave -glueto "$out"+tlrc # covariate mask
		    3dbucket masked_ISC_"$cov_n"_mask-"$mask"_main_t.nii.gz -glueto "$out"+tlrc # covariate z value thresholded
		    3dbucket masked_ISC_"$cov_n"_mask-"$mask"_main.nii.gz -glueto "$out"+tlrc # covariate beta thresholded
				    
		    # declare sub-brick to be t values
		    3drefit -substatpar 1 fitt  "$out"+tlrc

		    # rename sub-bricks
		    3drefit -sublabel 0 "main_""$cov_n" -sublabel 1 "main_z_""$cov_n" -sublabel 2 "mask_cov" "$out"+tlrc

		    # plot main effect
		    fig=fig_ISC_"$cov_n"_mask-"$mask"_effect-main
	        op=8
		    if [[ "$cov" == *"curBeta"* ]]; then
		        func_range=1
		    else
		        func_range=0.15
		    fi
	        
	        if [[ "$mask" == *"Caudate"* ]]; then
	            montx=3
	        else
	            montx=1
	        fi
	        
	        @chauffeur_afni -ulay $template -olay "$out"+tlrc -set_subbricks 0 0 0 -cbar Reds_and_Blues_Inv -func_range $func_range -prefix $fig -pbar_saveim "$fig"_pbar -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx $montx -monty 1 -box_focus_slices $mask_ave

	        
	        #@chauffeur_afni -ulay $template -olay out_ISC_"$cov_n"_mask-"$mask"_main.nii.gz -set_subbricks 0 0 0 -func_range $func_range -cbar Reds_and_Blues_Inv -prefix $fig -pbar_saveim "$fig"_pbar -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx $montx -monty 1 -box_focus_slices $mask_ave
        
        fi

    done

    # remove all ISC maps and clust stim output from directory
    rm mask* clust*
	   
done

#### at the end, create one plot showing all ISC in the ROIs

# add up all ISC maps in the ROIs
3dcalc -a out_ISC_"$task"_aHPC+tlrc.[0] -b out_ISC_"$task"_VTA+tlrc.[0] -c out_ISC_"$task"_NAcc+tlrc.[0] -d out_ISC_"$task"_Caudate+tlrc.[0] -prefix out_ISC_"$task"_mask-ROIs_main.nii.gz -expr 'a+b+c+d'

# create an image with all CMLE in ROIs
3dcalc -a out_ISC_fullCMLE_mask-aHPC+tlrc.[0] -b out_ISC_fullCMLE_mask-Caudate+tlrc.[0] -c out_ISC_fullCMLE_mask-NAcc+tlrc. -d out_ISC_fullCMLE_mask-VTA+tlrc.[0] -prefix out_ISC_fullCMLE_mask-ROIs_main.nii.gz -expr 'a+b+c+d'
3dcalc -a out_ISC_fullCMLE_mask-ROIs_main.nii.gz -expr 'isnegative(a)' -prefix mask_CMLE_ROI.nii.gz

3dcalc -a out_ISC_fullCMLE_mask-aHPC+tlrc.[0] -b out_ISC_fullCMLE_mask-NAcc+tlrc. -c out_ISC_fullCMLE_mask-VTA+tlrc.[0] -prefix out_ISC_fullCMLE_mask-ROIs_main.nii.gz -expr 'a+b+c+'
3dcalc -a out_ISC_fullCMLE_mask-ROIs_main.nii.gz -expr 'isnegative(a)' -prefix mask_CMLE_ROI.nii.gz

# create an image with all ROIs combined
3dcalc -a MNI_res-epi_label-aHPC_mask.nii.gz -b MNI_res-epi_label-VTASN_mask.nii.gz -c MNI_res-epi_label-NAcc_mask.nii.gz -d MNI_res-epi_label-Caudate_mask.nii.gz -prefix MNI_res-epi_label-ROIs_mask.nii.gz -expr 'a+b+c+d'

# plot ISC main effect
fig=fig_ISC_"$task"_mask-ROIs_effect-main
op=8
func_range=0.2
@chauffeur_afni -ulay $template -olay out_ISC_"$task"_mask-ROIs_main.nii.gz -set_subbricks 0 0 0 -func_range $func_range -cbar Reds_and_Blues_Inv -prefix $fig -pbar_saveim "$fig"_pbar -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY     

# plot CMLE main effect in ROIs
3dcalc -a out_ISC_fullCMLE_mask-aHPC+tlrc.[0] -b out_ISC_fullCMLE_mask-Caudate+tlrc.[0] -c out_ISC_fullCMLE_mask-NAcc+tlrc. -d out_ISC_fullCMLE_mask-VTA+tlrc.[0] -prefix out_ISC_fullCMLE_mask-ROIs_main.nii.gz -expr 'a+b+c+d'
3dcalc -a out_ISC_fullCMLE_mask-aHPC+tlrc.[2] -b out_ISC_fullCMLE_mask-NAcc+tlrc.[2] -c out_ISC_fullCMLE_mask-VTA+tlrc.[2] -prefix mask_CMLE_ROI.nii.gz -expr 'a+b+c' # do not include caudate to help with the slice positioning in the plot


fig=fig_ISC_CMLE_mask-ROIs_effect-main
op=8
func_range=1
@chauffeur_afni -ulay $template -olay out_ISC_fullCMLE_mask-ROIs_main.nii.gz -set_subbricks 0 0 0 -func_range $func_range -cbar Reds_and_Blues_Inv -prefix $fig -pbar_saveim "$fig"_pbar -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx 3 -monty 2 -box_focus_slices mask_CMLE_ROI.nii.gz     

# plot ROIs
fig=fig_ISC_"$task"_mask-ROIs
op=8
func_range=1
@chauffeur_afni -ulay $template -olay MNI_res-epi_label-ROIs_mask.nii.gz -set_subbricks 0 0 0 -func_range $func_range -cbar Reds_and_Blues_Inv -prefix $fig -pbar_saveim "$fig"_pbar -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 1 -box_focus_slices MNI_res-epi_label-ROIs_mask.nii.gz      


#### overlay main effects with Yeo networks

3dcalc -a $yeo -b out_ISC_curiosity_mask-GM+tlrc.[2] -expr 'a*b' -prefix yeo_curiosity.nii.gz
3dcalc -a $yeo -b out_ISC_memory_mask-GM+tlrc.[2] -expr 'a*b' -prefix yeo_memory.nii.gz

3dcalc -a $yeo -b out_ISC_fullCMLE_mask-GM+tlrc.[0] -expr 'ispositive(b)*a' -prefix yeo_clme_pos.nii.gz
3dcalc -a $yeo -b out_ISC_fullCMLE_mask-GM+tlrc.[0] -expr 'isnegative(b)*a' -prefix yeo_clme_neg.nii.gz

op=8
func_range=7
fig=fig_yeo_cmle_neg
@chauffeur_afni -ulay $template -olay yeo_clme_neg.nii.gz -set_subbricks 0 0 0 -func_range $func_range -cbar Color_circle_AJJ -prefix $fig -pbar_saveim "$fig"_pbar -pbar_posonly -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY       
fig=fig_yeo_cmle_pos
@chauffeur_afni -ulay $template -olay yeo_clme_pos.nii.gz -set_subbricks 0 0 0 -func_range $func_range -cbar Color_circle_AJJ -prefix $fig -pbar_saveim "$fig"_pbar -pbar_posonly -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY       

fig=fig_yeo_curiosity
@chauffeur_afni -ulay $template -olay yeo_curiosity.nii.gz -set_subbricks 0 0 0 -func_range $func_range -cbar Color_circle_AJJ -prefix $fig -pbar_saveim "$fig"_pbar -pbar_posonly -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY       

fig=fig_yeo_memory
@chauffeur_afni -ulay $template -olay yeo_memory.nii.gz -set_subbricks 0 0 0 -func_range $func_range -cbar Color_circle_AJJ -prefix $fig -pbar_saveim "$fig"_pbar -pbar_posonly -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY       

fig=fig_yeo
op=7
@chauffeur_afni -ulay $template -olay AFNI_Yeo2011_7Networks_MNI152_liberal.nii.gz -set_subbricks 0 0 0 -func_range $func_range -cbar Color_circle_AJJ -prefix $fig -pbar_saveim "$fig"_pbar -pbar_posonly -save_ftype "JPEG" -pbar_dim 64x1351H -opacity $op -label_mode 5 -label_color black -zerocolor white -montx 8 -monty 2 -box_focus_slices AMASK_FOCUS_ULAY    


