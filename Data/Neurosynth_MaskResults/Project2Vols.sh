#bash

for i in `seq 1 20`; do
wb_command -metric-to-volume-mapping IFD_NGR_Grads_Left/IFD777_NGR_Grad$(printf %02d $i).shape.gii ../../HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/S1200.L.midthickness_MSMAll.10k_fs_LR.surf.gii ../../HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/S1200_AverageT1w_restore.nii.gz IFD_NGR_Grads_Left/IFD777_NGR_Grad$(printf %02d $i)_Left.nii -ribbon-constrained ../../HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/S1200.L.white_MSMAll.10k_fs_LR.surf.gii ../../HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/S1200.L.pial_MSMAll.10k_fs_LR.surf.gii
done


for i in `seq 1 20`; do
wb_command -metric-to-volume-mapping IFD_NGR_Grads_Right/IFD777_NGR_Grad$(printf %02d $i).shape.gii ../../HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/S1200.R.midthickness_MSMAll.10k_fs_LR.surf.gii ../../HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/S1200_AverageT1w_restore.nii.gz IFD_NGR_Grads_Right/IFD777_NGR_Grad$(printf %02d $i)_Right.nii -ribbon-constrained ../../HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/S1200.R.white_MSMAll.10k_fs_LR.surf.gii ../../HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/S1200.R.pial_MSMAll.10k_fs_LR.surf.gii
done