Dists=load('/Volumes/eSSD/FinalCodes/CalcData/ALL777_IFD.mat');
IFD_NGR=Dists.IFdist_AVG{6,1};
IFD_GSR=Dists.IFdist_AVG{6,2};
IFD_AVG_NGR=mean(IFD_NGR);
IFD_AVG_GSR=mean(IFD_GSR);
clear Dists IFD_NGR IFD_GSR

% medial wall mask
mwL=gifti('/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/Human.MedialWall_Conte69.L.10k_fs_LR.label.gii');
mwR=gifti('/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/Human.MedialWall_Conte69.R.10k_fs_LR.label.gii');
HemiMask=logical((1-mwL.cdata') .* (1-mwR.cdata')); 

% surface for display
surfL=gifti('/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/S1200.L.inflated_MSMAll.10k_fs_LR.surf.gii');
avsurf_hemi.coord=surfL.vertices';
avsurf_hemi.tri=surfL.faces;

%% MetaLI

LItask=load('/Volumes/eSSD/Accessory/Data/metaLI_10k.mat');
LItask_m2=load('/Volumes/eSSD/Accessory/Data/metaLI_10k_method2.mat');
LItask_pos_md=LItask.metaposLIabs;
LItask_md=LItask.metaLIabs;

LItask2_md = LItask_m2.LImap; LItask2_md(isnan(LItask2_md))=0;
LItask2=gifti();
LItask2.cdata=LItask2_md';
save(LItask2,'/Volumes/eSSD/FinalCodes/ResultsMaps/LItask_method2.shape.gii','Base64Binary');

LI_10k=gifti('/Volumes/eSSD/Accessory/Data/LI_10k.shape.gii');
LI_10k_md=LI_10k.cdata(HemiMask);

SurfStatViewData_lxy( LItask_pos_md, avsurf_hemi, [0,0.8], '', 'Hori')
SurfStatViewData_lxy( LItask_md, avsurf_hemi, [0,0.8], '', 'Hori')
SurfStatViewData_lxy( LItask2_md, avsurf_hemi, [0,0.6], '', 'Hori')
SurfStatViewData_lxy( LItask_m2.termNum, avsurf_hemi, [0,295], '', 'Hori')

[prctile_LImeta_spin,r_LImeta_original,r_LImeta_rand]=spintest_10k(IFD_AVG_NGR',LItask2_md(HemiMask)',HemiMask,10000,'spin');
[prctile_LI10k_spin,r_LI10k_original,r_LI10k_rand]=spintest_10k(IFD_AVG_NGR',LI_10k_md',HemiMask,10000,'spin');


%% physical distance
% coordL=avsurf.coord(:,1:10242)';
% coordR=avsurf.coord(:,10243:end)';
% physicaldist=sqrt(sum((coordL - coordR) .^ 2,2))';
physd=gifti('/Volumes/eSSD/Accessory/Data/ToSurfIce/PhysicalDist_10k.shape.gii');
physd_md=physd.cdata(HemiMask);

[prctile_physd_spin,r_physd_original,r_physd_rand]=spintest_10k(IFD_AVG_NGR',physd_md,HemiMask,10000,'spin');

%% development
dev=gifti('/Volumes/eSSD/Accessory/Data/developmentAbsolute_10k.shape.gii');
dev_md=dev.cdata(HemiMask);

[prctile_dev_spin,r_dev_original,r_dev_rand]=spintest_10k(IFD_AVG_NGR',dev_md,HemiMask,10000,'spin');

%% evolution
evo=gifti('/Volumes/eSSD/Accessory/Data/evolutionAbsolute_10k.shape.gii');
evo_md=evo.cdata(HemiMask);

[prctile_evo_spin,r_evo_original,r_evo_rand]=spintest_10k(IFD_AVG_NGR',evo_md,HemiMask,10000,'spin');
% [prctile_evo_spin]=spinZtest_10k(FDGmapAVG(2,:)',evo_md,dev_md,HemiMask,10000,'spin');

%% myelin
myelin=gifti('/Volumes/eSSD/Accessory/Data/ToSurfIce/MyelinMap_S1200_10k.shape.gii');
myelin_md=myelin.cdata(HemiMask);

[prctile_myelin_spin,r_myelin_original,r_myelin_rand]=spintest_10k(IFD_AVG_NGR',myelin_md,HemiMask,10000,'spin');
