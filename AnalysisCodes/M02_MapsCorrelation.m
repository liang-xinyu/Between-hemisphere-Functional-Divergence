% Please change this path according to your saved path 
Codes_path = '/Volumes/eSSD/Between-hemisphere-Functional-Divergence/';

% This file need to be download separately 
Dists=load([Codes_path,'Data/ALL777_IFD.mat']);

IFD_NGR=Dists.IFdist_AVG{6,1};
IFD_AVG_NGR=mean(IFD_NGR);
clear Dists IFD_NGR IFD_GSR

% medial wall mask
mwL=gifti([Codes_path,'HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/Human.MedialWall_Conte69.L.10k_fs_LR.label.gii']);
mwR=gifti([Codes_path,'HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/Human.MedialWall_Conte69.R.10k_fs_LR.label.gii']);
HemiMask=logical((1-mwL.cdata') .* (1-mwR.cdata')); 

% surface for display
surfL=gifti([Codes_path,'HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/S1200.L.inflated_MSMAll.10k_fs_LR.surf.gii']);
avsurf_hemi.coord=surfL.vertices';
avsurf_hemi.tri=surfL.faces;

%% Spin test

% MetaLI
metaLI=gifti([Codes_path,'/HCP_Atlas/Corticalmaps/metaLI_tasks_10k.shape.gii']);
metaLI_md = metaLI.cdata(HemiMask);

[prctile_LImeta_spin,r_LImeta_original,r_LImeta_rand]=spintest_10k(IFD_AVG_NGR',metaLI_md,HemiMask,10000,'spin');

% development
dev=gifti([Codes_path,'/HCP_Atlas/Corticalmaps/developmentAbsolute_10k.shape.gii']);
dev_md=dev.cdata(HemiMask);

[prctile_dev_spin,r_dev_original,r_dev_rand]=spintest_10k(IFD_AVG_NGR',dev_md,HemiMask,10000,'spin');

% evolution
evo=gifti([Codes_path,'/HCP_Atlas/Corticalmaps/evolutionAbsolute_10k.shape.gii']);
evo_md=evo.cdata(HemiMask);

[prctile_evo_spin,r_evo_original,r_evo_rand]=spintest_10k(IFD_AVG_NGR',evo_md,HemiMask,10000,'spin');
% [prctile_evo_spin]=spinZtest_10k(FDGmapAVG(2,:)',evo_md,dev_md,HemiMask,10000,'spin');

% myelin
myelin=gifti([Codes_path,'/HCP_Atlas/Corticalmaps/MyelinMap_S1200_10k.shape.gii']);
myelin_md=myelin.cdata(HemiMask);

[prctile_myelin_spin,r_myelin_original,r_myelin_rand]=spintest_10k(IFD_AVG_NGR',myelin_md,HemiMask,10000,'spin');
