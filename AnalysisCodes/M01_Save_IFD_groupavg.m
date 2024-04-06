IFD=load('/Volumes/eSSD/FinalCodes/CalcData/ALL777_IFD.mat');
% IFD in 6D
IFdist6_NGR=IFD.IFdist_AVG{6,1}; %NGR
IFdist6_GSR=IFD.IFdist_AVG{6,2}; %GSR
IFdist6_ICC_NGR=IFD.ICC_IFdist{6,1}; %NGR
IFdist6_ICC_GSR=IFD.ICC_IFdist{6,2}; %GSR


mwL=gifti('/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/Human.MedialWall_Conte69.L.10k_fs_LR.label.gii');
mwR=gifti('/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/Human.MedialWall_Conte69.R.10k_fs_LR.label.gii');
HemiMask=logical((1-mwL.cdata') .* (1-mwR.cdata')); 

% surface for display
OutPath = '/Volumes/eSSD/FinalCodes/CalcData/SurfMetrics';

IFDavg_NGR=gifti();
Tcdata = zeros(10242,1);Tcdata(HemiMask)=mean(IFdist6_NGR)';
IFDavg_NGR.cdata=Tcdata;
save(IFDavg_NGR,[OutPath,filesep,'IFD_Avg_NGR_ALL777.shape.gii'],'Base64Binary');

IFDavg_GSR=gifti();
Tcdata = zeros(10242,1);Tcdata(HemiMask)=mean(IFdist6_GSR)';
IFDavg_GSR.cdata=Tcdata;
save(IFDavg_GSR,[OutPath,filesep,'IFD_Avg_GSR_ALL777.shape.gii'],'Base64Binary');

IFD_ICC_NGR=gifti();
Tcdata = zeros(10242,1);Tcdata(HemiMask)=IFdist6_ICC_NGR';
IFD_ICC_NGR.cdata=Tcdata;
save(IFD_ICC_NGR,[OutPath,filesep,'IFD_ICC_NGR_ALL777.shape.gii'],'Base64Binary');

IFD_ICC_GSR=gifti();
Tcdata = zeros(10242,1);Tcdata(HemiMask)=IFdist6_ICC_GSR';
IFD_ICC_GSR.cdata=Tcdata;
save(IFD_ICC_GSR,[OutPath,filesep,'IFD_ICC_GSR_ALL777.shape.gii'],'Base64Binary');

