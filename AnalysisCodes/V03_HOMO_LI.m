mwL=gifti('/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/Human.MedialWall_Conte69.L.10k_fs_LR.label.gii');
mwR=gifti('/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/Human.MedialWall_Conte69.R.10k_fs_LR.label.gii');
HemiMask=logical((1-mwL.cdata') .* (1-mwR.cdata')); 

% surface for display
surfL=gifti('/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/S1200.L.inflated_MSMAll.10k_fs_LR.surf.gii');
avsurf_hemi.coord=surfL.vertices';
avsurf_hemi.tri=surfL.faces;

IDfile=load('/Volumes/eSSD/FinalCodes/Inform/Subinfo_rhand_nocode.mat');
conID=IDfile.Sub_ID;

R1_NGR_path=g_ls('/Volumes/eSSD/Accessory/Data/REST1/Vertex_Homo/*_REST1_NGR_hemi_FC_Z_HomoxFCS.mat');
R2_NGR_path=g_ls('/Volumes/eSSD/Accessory/Data/REST2/Vertex_Homo/*_REST2_NGR_hemi_FC_Z_HomoxFCS.mat');


si1=cellfun(@(x) str2double(x(end-38:end-33)),R1_NGR_path);
si2=cellfun(@(x) str2double(x(end-38:end-33)),R2_NGR_path);
[~,is1a,is1b] = intersect(conID,si1);
[~,is2a,is2b] = intersect(conID,si2);

Homo_R1=zeros(777,9354);
Homo_R2=zeros(777,9354);
LIiFCS_R1=zeros(777,9354);
LIiFCS_R2=zeros(777,9354);
for i=1:length(conID)
    tmpg1=load(R1_NGR_path{is1b(i)});
    tmpg2=load(R2_NGR_path{is2b(i)});
    Homo_R1(i,:)=tmpg1.Homotopy(HemiMask);
    Homo_R2(i,:)=tmpg2.Homotopy(HemiMask);
    LIiFCS_R1(i,:)=abs(tmpg1.FCS_pos_L(HemiMask)-tmpg1.FCS_pos_R(HemiMask))./(tmpg1.FCS_pos_L(HemiMask)+tmpg1.FCS_pos_R(HemiMask));
    LIiFCS_R2(i,:)=abs(tmpg2.FCS_pos_L(HemiMask)-tmpg2.FCS_pos_R(HemiMask))./(tmpg2.FCS_pos_L(HemiMask)+tmpg2.FCS_pos_R(HemiMask));
end

Homo_global_R1=mean(Homo_R1,2);
Homo_global_R2=mean(Homo_R2,2);
LIiFCS_global_R1=mean(LIiFCS_R1,2);
LIiFCS_global_R2=mean(LIiFCS_R2,2);

ICC_Homo_whole=IPN_icc([Homo_global_R1,Homo_global_R2],2,'single');
ICC_LIiFCS_whole=IPN_icc([LIiFCS_global_R1,LIiFCS_global_R2],2,'single');
for j = 1:9354
    ICC_Homo_vertex(1,j) = IPN_icc([Homo_R1(:,j),Homo_R2(:,j)],2,'single');
    ICC_LIiFCS_vertex(1,j) = IPN_icc([LIiFCS_R1(:,j),LIiFCS_R2(:,j)],2,'single');
end

HoFC_groupavg=zeros(1,10242);
HoFC_groupavg(HemiMask)=mean((Homo_R1+Homo_R2)/2);
HoFCavg=gifti();
HoFCavg.cdata=HoFC_groupavg';
save(HoFCavg,'/Volumes/eSSD/FinalCodes/ResultsMaps/HoFC777_NGR_AVG.shape.gii','Base64Binary');

LIiFCSFC_groupavg=zeros(1,10242);
LIiFCSFC_groupavg(HemiMask)=mean((LIiFCS_R1+LIiFCS_R2)/2);
LIiFCSavg=gifti();
LIiFCSavg.cdata=LIiFCSFC_groupavg';
save(LIiFCSavg,'/Volumes/eSSD/FinalCodes/ResultsMaps/LIiFCS777_NGR_AVG.shape.gii','Base64Binary');

HoFC_ICC=zeros(1,10242);
HoFC_ICC(HemiMask)=ICC_Homo_vertex;
HoFCICC=gifti();
HoFCICC.cdata=HoFC_ICC';
save(HoFCICC,'/Volumes/eSSD/FinalCodes/ResultsMaps/HoFC777_NGR_ICC.shape.gii','Base64Binary');

LIiFCSFC_ICC=zeros(1,10242);
LIiFCSFC_ICC(HemiMask)=ICC_LIiFCS_vertex;
LIiFCSICC=gifti();
LIiFCSICC.cdata=LIiFCSFC_ICC';
save(LIiFCSICC,'/Volumes/eSSD/FinalCodes/ResultsMaps/LIiFCS777_NGR_ICC.shape.gii','Base64Binary');


mean(ICC_Homo_vertex),std(ICC_Homo_vertex)
mean(ICC_LIiFCS_vertex),std(ICC_LIiFCS_vertex)
%%
HeadMotion = load('/Volumes/eSSD/Project-HCP-HemiSurfAsym/GradientAna/HeadMotion_FD_777subs.mat');
[~,~,ihm]=intersect(conID,HeadMotion.conID,'stable');
mFD = (HeadMotion.mFD1(ihm)+HeadMotion.mFD2(ihm))/2;
%mFD_scrub = (HeadMotion.mFD1_scrub(ihm)+HeadMotion.mFD2_scrub(ihm))/2;

% cognition
[numericData, textData, ~]=xlsread('/Volumes/eSSD/Accessory/Homotopy-IntraLI-Pipeline/Utilized/HCP_Cognition_byLXY0318.xlsx');
[~,~,cogia]=intersect(conID,numericData(:,1),'stable');
[Dem_numData, Dem_textData, ~]=xlsread('/Volumes/eSSD/Accessory/Data/Demographics_for_IQ_0714.xlsx');
[~,~,demia]=intersect(conID,Dem_numData(:,1),'stable');

%% demographic
DAge=Dem_numData(demia,2);
DSex=Dem_textData(demia+1,7);
DSexN=strcmp(DSex,'F')+1; % num type, 1=male, 2=female
DHand=Dem_numData(demia,5);
DICV=zscore(Dem_numData(demia,8));

%% FIQ scores from NIH

FIQ_cols = [48,4,6,8,17,46];
Gf_cogs=numericData(cogia,FIQ_cols);
Gf_cogs_titles=textData(1,FIQ_cols)';
Vfd=~isnan(sum(Gf_cogs,2));
sum(Vfd)

Homo_global=(Homo_global_R1+Homo_global_R2)/2;
LIiFCS_global=(LIiFCS_global_R1+LIiFCS_global_R2)/2;

Cogouts = isoutlier(Gf_cogs(Vfd,:),"mean");
Allids=1:777;
VFids=Allids(Vfd);
Vfd(VFids(logical(sum(Cogouts,2))))=0;
sum(Vfd)

IFDouts=logical(sum(isoutlier([Homo_global,LIiFCS_global],"mean"),2));
Vfd(IFDouts)=0;
sum(Vfd)

ICVouts=isoutlier(DICV,"mean");
Vfd(ICVouts)=0;
sum(Vfd)

[r_IFDxICV(1),p_IFDxICV(1)] = partialcorr(Homo_global(Vfd),DICV(Vfd),[DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
[r_IFDxICV(2),p_IFDxICV(2)] = partialcorr(LIiFCS_global(Vfd),DICV(Vfd),[DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
for i=1:6
    [rw(i,1),pw(i,1)] = partialcorr(Homo_global(Vfd),Gf_cogs(Vfd,i),[DICV(Vfd),DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
    [rw(i,2),pw(i,2)] = partialcorr(LIiFCS_global(Vfd),Gf_cogs(Vfd,i),[DICV(Vfd),DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
end
