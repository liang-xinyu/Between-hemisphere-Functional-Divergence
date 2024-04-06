R1=load('/Volumes/eSSD/Accessory/Data/GroupAVG/Template_777_Com/GroupAVG_777_R1NGR_gradient.mat');
R2=load('/Volumes/eSSD/Accessory/Data/GroupAVG/Template_777_Com/GroupAVG_777_R2NGR_gradient.mat');

explain_R1=R1.g_hemi.lambda{1}(1:10)/sum(R1.g_hemi.lambda{1}(1:10));
explain_R2=R2.g_hemi.lambda{1}(1:10)/sum(R2.g_hemi.lambda{1}(1:10));

explain_R1_all=R1.g_hemi.lambda{1}/sum(R1.g_hemi.lambda{1});
explain_R2_all=R2.g_hemi.lambda{1}/sum(R2.g_hemi.lambda{1});


Dists=load('/Volumes/eSSD/FinalCodes/CalcData/ALL777_IFD.mat');
ICC=Dists.ICC_IFdist_whole(:,1);

explain10=(explain_R1+explain_R2)/2;
explainAlL=(explain_R1_all+explain_R2_all)/2;

mwL=gifti('/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/Human.MedialWall_Conte69.L.10k_fs_LR.label.gii');
mwR=gifti('/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/Human.MedialWall_Conte69.R.10k_fs_LR.label.gii');
HemiMask=logical((1-mwL.cdata') .* (1-mwR.cdata')); 

% surface for display
OutPath = '/Volumes/eSSD/FinalCodes/ResultsMaps/';

IFDavg_NGR=gifti();
Tcdata = zeros(10242,1);Tcdata(HemiMask)=Dists.ICC_IFdist{6,1}(1,:)';
IFDavg_NGR.cdata=Tcdata;
save(IFDavg_NGR,[OutPath,filesep,'IFD_Sub1_NGR_R1.shape.gii'],'Base64Binary');