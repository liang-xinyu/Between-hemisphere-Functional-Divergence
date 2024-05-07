% Please change this path according to your saved path 
Codes_path = '/Volumes/eSSD/Between-hemisphere-Functional-Divergence/';

% load data
% This file need to be download separately 
Dists=load([Codes_path,'Data/ALL777_IFD.mat']);

IFDw6_NGR=Dists.IFdistAll_AVG{1}(:,6); %NGR
IFDw6_GSR=Dists.IFdistAll_AVG{2}(:,6); %GSR

conID=Dists.SubID;
clear Dists

% medial wall mask
mwL=gifti([Codes_path,'HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/Human.MedialWall_Conte69.L.10k_fs_LR.label.gii']);
mwR=gifti([Codes_path,'HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/Human.MedialWall_Conte69.R.10k_fs_LR.label.gii']);
HemiMask=logical((1-mwL.cdata') .* (1-mwR.cdata')); 

% surface for display
surfL=gifti([Codes_path,'HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/S1200.L.inflated_MSMAll.10k_fs_LR.surf.gii']);
avsurf_hemi.coord=surfL.vertices';
avsurf_hemi.tri=surfL.faces;

%% Varibles
% subject demographics
SubSize=length(conID);

HeadMotion = load([Codes_path,'Data/HeadMotion_FD_777subs.mat']);
[~,~,ihm]=intersect(conID,HeadMotion.conID,'stable');
mFD = (HeadMotion.mFD1(ihm)+HeadMotion.mFD2(ihm))/2;

% cognition
[numericData, textData, ~]=xlsread([Codes_path,'Data/HCP_Cognition_byLXY0318.xlsx']);
[~,~,cogia]=intersect(conID,numericData(:,1),'stable');
[Dem_numData, Dem_textData, ~]=xlsread([Codes_path,'Data/Demographics_for_IQ_0714.xlsx']);
[~,~,demia]=intersect(conID,Dem_numData(:,1),'stable');

% demographic
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

% outliers
Cogouts = isoutlier(Gf_cogs(Vfd,:),"mean");
Allids=1:777;
VFids=Allids(Vfd);
Vfd(VFids(logical(sum(Cogouts,2))))=0;
sum(Vfd)

IFDouts=logical(sum(isoutlier([IFDw6_NGR,IFDw6_GSR],"mean"),2));
Vfd(IFDouts)=0;
sum(Vfd)

ICVouts=isoutlier(DICV,"mean");
Vfd(ICVouts)=0;
sum(Vfd)

% relationship between IFD and ICV
[r_IFDxICV(1),p_IFDxICV(1)] = partialcorr(IFDw6_NGR(Vfd),DICV(Vfd),[DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
[r_IFDxICV(2),p_IFDxICV(2)] = partialcorr(IFDw6_GSR(Vfd),DICV(Vfd),[DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);

% relationship between IFD and Age
[r,p] = partialcorr(IFDw6_NGR(Vfd),DAge(Vfd),[DICV(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
[r,p] = partialcorr(IFDw6_GSR(Vfd),DAge(Vfd),[DICV(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);

% differences of IFD between Sex
SurfStatPlot(DSex(Vfd),IFDw6_NGR(Vfd),1+term(DAge(Vfd))+term(DHand(Vfd))+term(DICV(Vfd))+term(mFD(Vfd)));
SurfStatPlot(DSex(Vfd),IFDw6_GSR(Vfd),1+term(DAge(Vfd))+term(DHand(Vfd))+term(DICV(Vfd))+term(mFD(Vfd)));

% relationship between IFD and Handedness
[r,p] = partialcorr(IFDw6_NGR(Vfd),DHand(Vfd),[DICV(Vfd),DAge(Vfd),DSexN(Vfd),mFD(Vfd)]);
[r,p] = partialcorr(IFDw6_GSR(Vfd),DHand(Vfd),[DICV(Vfd),DAge(Vfd),DSexN(Vfd),mFD(Vfd)]);

% relationship between IFD and cognitions
for i=1:6
    [rw(i,1),pw(i,1)] = partialcorr(IFDw6_NGR(Vfd),Gf_cogs(Vfd,i),[DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
    [rw(i,2),pw(i,2)] = partialcorr(IFDw6_GSR(Vfd),Gf_cogs(Vfd,i),[DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);    
    [rw(i,3),pw(i,3)] = partialcorr(IFDw6_NGR(Vfd),Gf_cogs(Vfd,i),[DICV(Vfd),DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
    [rw(i,4),pw(i,4)] = partialcorr(IFDw6_GSR(Vfd),Gf_cogs(Vfd,i),[DICV(Vfd),DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
    [r_ICVxCog(i),p_ICVxCog(i)] = partialcorr(DICV(Vfd),Gf_cogs(Vfd,i),[DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
end

%% global mediation
% need mediation_toolbox
X = zscore(DICV(Vfd));
Y = zscore(Gf_cogs(Vfd,1));
M = zscore(IFDw6_NGR(Vfd));
covs = zscore([DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
[paths, stats] = mediation(X, Y, M, 'covs',covs, 'boot', 'verbose', 'bootsamples', 10000, 'doCIs');

X = zscore(DICV(Vfd));
Y = zscore(Gf_cogs(Vfd,3));
M = zscore(IFDw6_NGR(Vfd));
covs = zscore([DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
[paths, stats] = mediation(X, Y, M, 'covs',covs, 'boot', 'verbose', 'bootsamples', 10000, 'doCIs');

X = zscore(DICV(Vfd));
Y = zscore(Gf_cogs(Vfd,1));
M = zscore(IFDw6_GSR(Vfd));
covs = zscore([DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
[paths, stats] = mediation(X, Y, M, 'covs',covs, 'boot', 'verbose', 'bootsamples', 10000, 'doCIs');

X = zscore(DICV(Vfd));
Y = zscore(Gf_cogs(Vfd,3));
M = zscore(IFDw6_GSR(Vfd));
covs = zscore([DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
[paths, stats] = mediation(X, Y, M, 'covs',covs, 'boot', 'verbose', 'bootsamples', 10000, 'doCIs');

