#!/bin/bash

################################################################################
# Check WM and Vent masks and ROI masks
#https://neurovault.org/collections/12980/
#Save ROIs: 
#	MNI_res-epi_label-aHPC_mask.nii.gz 
#	MNI_res-epi_label-Caudate_mask.nii.gz 
#	MNI_res-epi_label-NAcc_mask.nii.gz 
#	MNI_res-epi_label-VTASN_mask.nii.gz
#	in 'outdir=$derivroot/ROI' directory
################################################################################

# Set top level directory structure
cd ~
topdir=/mnt/dell_storage/labs/atsuchiyagaito/MIL
echo $topdir

derivroot=$topdir/derivatives
outdir=$derivroot/ROI
mkdir -p $outdir

# define subject listecho $
BIDSdir=$topdir/rawdata

cd $BIDSdir


# define tasks to loop through
tasks=(magictrickwatching rest)
# define subjects to loop through
#subjects=(`ls -d sub*`)
subjects=("sub-control001" "sub-control002" "sub-control003")


#Combine ROIs
3dcalc -a ${outdir}/MNI_res-epi_label-aHPC_mask.nii.gz -b ${outdir}/MNI_res-epi_label-Caudate_mask.nii.gz \
	-c ${outdir}/MNI_res-epi_label-NAcc_mask.nii.gz -d ${outdir}/MNI_res-epi_label-VTASN_mask.nii.gz \
	-expr 'a + b + c + d' -prefix ${outdir}/rm.MNI_res-epi_label_ROIs_mask.nii.gz
3dcalc -a ${outdir}/rm.MNI_res-epi_label_ROIs_mask.nii.gz -expr 'equals(a,1)' -prefix ${outdir}/MNI_res-epi_label_ROIs_mask.nii.gz


# for each task
for task in "${tasks[@]}"; do

	# initialize an empty array to hold the subject filenames for this task
	subj_files=()

	# for each subject in the subjects array
	for subj in "${subjects[@]}"; do
		# determine subjstr
		subjstr="${subj}_task-${task}"

		# add this subject's filename to the array
		subj_files_gm+=("${derivroot}/${subj}/func/${subjstr}_space-MNI152NLin2009cAsym_label-epiGM_mask.nii.gz")
		subj_files_wm+=("${derivroot}/${subj}/func/${subjstr}_space-MNI152NLin2009cAsym_label-epiWM_mask.nii.gz")
		subj_files_vent+=("${derivroot}/${subj}/func/${subjstr}_space-MNI152NLin2009cAsym_label-epiVent_mask.nii.gz")
	done

	#GM
	# now call 3dMean with all the subject filenames for this task
	3dMean -datum float -prefix ${outdir}/${task}_gm_average.nii.gz "${subj_files_gm[@]}"
	#creating mask
	3dcalc -datum float -prefix ${outdir}/${task}_gm_mask.nii.gz -a ${outdir}/${task}_gm_average.nii.gz -expr 'ispositive(a-0.499)'
	
	#WM
	# now call 3dMean with all the subject filenames for this task
	3dMean -datum float -prefix ${outdir}/${task}_wm_average.nii.gz "${subj_files_wm[@]}"
	#creating mask
	3dcalc -datum float -prefix ${outdir}/${task}_wm_mask.nii.gz -a ${outdir}/${task}_wm_average.nii.gz -expr 'ispositive(a-0.499)'
	
	#Vent
	# now call 3dMean with all the subject filenames for this task
	3dMean -datum float -prefix ${outdir}/${task}_vent_average.nii.gz "${subj_files_vent[@]}"
	#creating mask
	3dcalc -datum float -prefix ${outdir}/${task}_vent_mask.nii.gz -a ${outdir}/${task}_vent_average.nii.gz -expr 'ispositive(a-0.499)'
	
	#Combine all
	#3dcalc -a ${outdir}/${task}_gm_mask.nii.gz -b ${outdir}/${task}_wm_mask.nii.gz -c ${outdir}/${task}_vent_mask.nii.gz -expr '1*a + 2*b + 4*c' -prefix ${outdir}/${task}_gm_wm_vent_mask.nii.gz
	
	#Combine WM and Vent
	3dcalc -a ${outdir}/${task}_wm_mask.nii.gz -b ${outdir}/${task}_vent_mask.nii.gz -expr 'a + b ' -prefix ${outdir}/rm.${task}_wm_vent_mask.nii.gz
	3dcalc -a ${outdir}/rm.${task}_wm_vent_mask.nii.gz -expr 'equals(a,1)' -prefix ${outdir}/${task}_wm_vent_mask.nii.gz
	
	#Create figures
	fig=${outdir}/fig_${task}_ROIs_wm_vent
	op=5
	ulay=${outdir}/${task}_wm_vent_mask.nii.gz
	olay=${outdir}/MNI_res-epi_label_ROIs_mask.nii.gz
	@chauffeur_afni -ulay $ulay -olay $olay \
		-set_subbricks 0 0 0 -prefix $fig -save_ftype "JPEG" \
		-opacity $op -cbar "red_monochrome" -label_color white -zerocolor black \
		-montx 8 -monty 1 -box_focus_slices $olay -label_size 4
done

rm ${outdir}/rm*






