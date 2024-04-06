% IFD dim6
Dists=load('/Volumes/eSSD/FinalCodes/CalcData/ALL777_IFD.mat');
IFD_NGR=Dists.IFdist_AVG{6,1};
IFD_GSR=Dists.IFdist_AVG{6,2};

% medial wall mask
mwL=gifti('/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/Human.MedialWall_Conte69.L.10k_fs_LR.label.gii');
mwR=gifti('/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/Human.MedialWall_Conte69.R.10k_fs_LR.label.gii');
HemiMask=logical((1-mwL.cdata') .* (1-mwR.cdata')); 
MW=gifti();
MW.cdata=~HemiMask';
save(MW,'/Volumes/eSSD/FinalCodes/ResultsMaps/MedialWall_LRcombined_10k.shape.gii','Base64Binary');


% surface for display
surfL=gifti('/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/S1200.L.inflated_MSMAll.10k_fs_LR.surf.gii');
surfR=gifti('/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/S1200.R.inflated_MSMAll.10k_fs_LR.surf.gii');
avsurf.coord=[surfL.vertices;surfR.vertices]';
avsurf.tri=[surfL.faces;surfR.faces+10242];
avsurf_hemi.coord=surfL.vertices';
avsurf_hemi.tri=surfL.faces;

IFD_ICC=zeros(1,10242);
IFD_ICC(HemiMask)=Dists.ICC_IFdist{6,1};
IFDICC=gifti();
IFDICC.cdata=IFD_ICC';
save(IFDICC,'/Volumes/eSSD/FinalCodes/ResultsMaps/IFD777_NGR_ICC.shape.gii','Base64Binary');

IFD_GSR_ICC=zeros(1,10242);
IFD_GSR_ICC(HemiMask)=Dists.ICC_IFdist{6,2};
IFDICC_GSR=gifti();
IFDICC_GSR.cdata=IFD_GSR_ICC';
save(IFDICC_GSR,'/Volumes/eSSD/FinalCodes/ResultsMaps/IFD777_GSR_ICC.shape.gii','Base64Binary');



IFD_groupavg=zeros(1,10242);
IFD_groupavg(HemiMask)=mean(IFD_NGR);
IFD_groupstd=zeros(1,10242);
IFD_groupstd(HemiMask)=std(IFD_NGR);
SurfStatViewData_lxy( IFD_groupavg, avsurf_hemi, [0,0.1], '', 'Hori')
SurfStatViewData_lxy( IFD_groupstd, avsurf_hemi, [0,0.05], '', 'Hori')
IFD_CVs=zeros(1,10242);
IFD_CVs(HemiMask)=std(IFD_NGR)./mean(IFD_NGR);
SurfStatViewData_lxy( IFD_CVs, avsurf_hemi, [0,0.8], '', 'Hori')

IFDavg=gifti();
IFDavg.cdata=IFD_groupavg';
save(IFDavg,'/Volumes/eSSD/FinalCodes/ResultsMaps/IFD777_NGR_AVG.shape.gii','Base64Binary');
IFDstd=gifti();
IFDstd.cdata=IFD_groupstd';
save(IFDstd,'/Volumes/eSSD/FinalCodes/ResultsMaps/IFD777_NGR_SD.shape.gii','Base64Binary');

IFD_groupavg_GSR=zeros(1,10242);
IFD_groupavg_GSR(HemiMask)=mean(IFD_GSR);
IFDavg_GSR=gifti();
IFDavg_GSR.cdata=IFD_groupavg_GSR';
save(IFDavg_GSR,'/Volumes/eSSD/FinalCodes/ResultsMaps/IFD777_GSR_AVG.shape.gii','Base64Binary');


Pavg = prctile(mean(IFD_NGR),5:5:95);
IFD_grades = discretize(mean(IFD_NGR),[0,Pavg,1]);
IFD_grades_map=zeros(20,10242);
for i=1:20
    IFD_grades_map(i,HemiMask)=IFD_grades==i;
    IFD_gradmaps=gifti();
    IFD_gradmaps.cdata=IFD_grades_map(i,:)';
    save(IFD_gradmaps,['/Volumes/eSSD/FinalCodes/ResultsMaps/IFD777_NGR_Grad',num2str(i,'%02d'),'.shape.gii'],'Base64Binary');
end
SurfStatViewData_lxy( IFD_grades_map(20,:), avsurf_hemi, [0,1], '', 'Hori')

