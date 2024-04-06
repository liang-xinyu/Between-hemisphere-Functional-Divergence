R1_NGR_Grads=load('/Volumes/eSSD/Accessory/Data/GroupAVG/Template_777_Com/GroupAVG_777_R1NGR_gradient.mat');
R2_NGR_Grads=load('/Volumes/eSSD/Accessory/Data/GroupAVG/Template_777_Com/GroupAVG_777_R2NGR_gradient.mat');

% medial wall mask
mwL=gifti('/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/Human.MedialWall_Conte69.L.10k_fs_LR.label.gii');
mwR=gifti('/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/Human.MedialWall_Conte69.R.10k_fs_LR.label.gii');
HemiMask=logical((1-mwL.cdata') .* (1-mwR.cdata')); 

% surface for display
surfL=gifti('/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/S1200.L.inflated_MSMAll.10k_fs_LR.surf.gii');
avsurf_hemi.coord=surfL.vertices';
avsurf_hemi.tri=surfL.faces;

for i=1:10
    GiftiTemp1=gifti();
    GTdata1=zeros(10242,1);
    GTdata1(HemiMask)=R1_NGR_Grads.g_hemi.gradients{1}(:,i);
    GiftiTemp1.cdata=GTdata1;
    save(GiftiTemp1,['/Volumes/eSSD/FinalCodes/ResultsMaps/TemplateGradients/ALL777_Gradient',num2str(i,'%2d'),'_map_R1_NGR.shape.gii'],'Base64Binary');
    figure,SurfStatViewData_lxy( GTdata1*HemiMask, avsurf_hemi, [min(GTdata1),max(GTdata1)], '', 'Hori')

    GiftiTemp2=gifti();
    GTdata2=zeros(10242,1);
    GTdata2(HemiMask)=R2_NGR_Grads.g_hemi.gradients{1}(:,i);
    GiftiTemp2.cdata=GTdata2;
    save(GiftiTemp2,['/Volumes/eSSD/FinalCodes/ResultsMaps/TemplateGradients/ALL777_Gradient',num2str(i,'%2d'),'_map_R2_NGR.shape.gii'],'Base64Binary');
    figure,SurfStatViewData_lxy( GTdata2*HemiMask, avsurf_hemi, [min(GTdata2),max(GTdata2)], '', 'Hori')

end

