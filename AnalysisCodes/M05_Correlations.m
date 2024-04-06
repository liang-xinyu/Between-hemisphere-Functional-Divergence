% load data (Distance in diverse dimensions)
Dists=load('/Volumes/eSSD/FinalCodes/CalcData/ALL777_IFD.mat');
% Dists=load('/Volumes/eSSD/FinalCodes/CalcData/ALL777_R2L_IFD.mat');

% Dists=load('/Volumes/eSSD/FinalCodes/CalcData/Uri387_IFD.mat');
% IFD in 6D

IFDw6_NGR=Dists.IFdistAll_AVG{1}(:,6); %NGR
IFDw6_GSR=Dists.IFdistAll_AVG{2}(:,6); %GSR
% IFD_NGR=Dists.IFdist_AVG{6,1}; %NGR
% IFD_GSR=Dists.IFdist_AVG{6,2}; %GSR

conID=Dists.SubID;
clear Dists

% medial wall mask
mwL=gifti('/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/Human.MedialWall_Conte69.L.10k_fs_LR.label.gii');
mwR=gifti('/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/Human.MedialWall_Conte69.R.10k_fs_LR.label.gii');
HemiMask=logical((1-mwL.cdata') .* (1-mwR.cdata')); 

% surface for display
surfL=gifti('/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/S1200.L.inflated_MSMAll.10k_fs_LR.surf.gii');
avsurf_hemi.coord=surfL.vertices';
avsurf_hemi.tri=surfL.faces;

%% smooth
% [IFD_NGR_s8_name] = SurfSmooth_wb(IFD_NGR,8,HemiMask,'IFD_NGR_6Dinds_Uri387');
% [IFD_GSR_s8_name] = SurfSmooth_wb(IFD_GSR,8,HemiMask,'IFD_GSR_6Dinds_Uri387');
% IFD_NGR_s8_data=gifti(IFD_NGR_s8_name);
% IFD_GSR_s8_data=gifti(IFD_GSR_s8_name);

%% load data
% surfdata = zeros(777,10242);
% surfdata(:,HemiMask) = IFD_NGR;
% surfdata_GSR = zeros(777,10242);
% surfdata_GSR(:,HemiMask) = IFD_GSR;
load('/Volumes/eSSD/FinalCodes/IFD_6Dinds_s8.mat');

%% varis
% subject demographics
SubSize=length(conID);
% Subinfo=load('/Volumes/eSSD/FinalCodes/Inform/Subinfo_rhand_nocode.mat');
% SubID=Subinfo.Sub_ID;
% [~,~,ics]=intersect(conID,SubID);
HeadMotion = load('/Volumes/eSSD/Project-HCP-HemiSurfAsym/GradientAna/HeadMotion_FD_777subs.mat');
[~,~,ihm]=intersect(conID,HeadMotion.conID,'stable');
mFD = (HeadMotion.mFD1(ihm)+HeadMotion.mFD2(ihm))/2;
%mFD_scrub = (HeadMotion.mFD1_scrub(ihm)+HeadMotion.mFD2_scrub(ihm))/2;

%% cognition
[numericData, textData, ~]=xlsread('/Volumes/eSSD/Accessory/Homotopy-IntraLI-Pipeline/Utilized/HCP_Cognition_byLXY0318.xlsx');
[~,~,cogia]=intersect(conID,numericData(:,1),'stable');
[Dem_numData, Dem_textData, ~]=xlsread('/Volumes/eSSD/Accessory/Data/Demographics_for_IQ_0714.xlsx');
[~,~,demia]=intersect(conID,Dem_numData(:,1),'stable');

%% demographic
DAge=Dem_numData(demia,2);
DSex=Dem_textData(demia+1,7);
DSexN=strcmp(DSex,'F')+1; % num type, 1=male, 2=female
DHand=Dem_numData(demia,5);
%DEduc=Dem_numData(demia,6);
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

[r_IFDxICV(1),p_IFDxICV(1)] = partialcorr(IFDw6_NGR(Vfd),DICV(Vfd),[DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
[r_IFDxICV(2),p_IFDxICV(2)] = partialcorr(IFDw6_GSR(Vfd),DICV(Vfd),[DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);

[r,p] = partialcorr(IFDw6_NGR(Vfd),DAge(Vfd),[DICV(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
[r,p] = partialcorr(IFDw6_GSR(Vfd),DAge(Vfd),[DICV(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);

[ t, df, pval, SSS ] = SurfStatPlot(DSex(Vfd),IFDw6_NGR(Vfd),1+term(DAge(Vfd))+term(DHand(Vfd))+term(DICV(Vfd))+term(mFD(Vfd)))
SurfStatPlot(DSex(Vfd),IFDw6_GSR(Vfd),1+term(DAge(Vfd))+term(DHand(Vfd))+term(DICV(Vfd))+term(mFD(Vfd)))
x=DSex(Vfd);
save('SexxIFD6_NGR.mat','x','SSS');

[ t, df, pval, SSS ] = SurfStatPlot(DSex(Vfd),IFDw6_GSR(Vfd),1+term(DAge(Vfd))+term(DHand(Vfd))+term(DICV(Vfd))+term(mFD(Vfd)))

[r,p] = partialcorr(IFDw6_NGR(Vfd),DHand(Vfd),[DICV(Vfd),DAge(Vfd),DSexN(Vfd),mFD(Vfd)]);
[r,p] = partialcorr(IFDw6_GSR(Vfd),DHand(Vfd),[DICV(Vfd),DAge(Vfd),DSexN(Vfd),mFD(Vfd)]);

for i=1:6
    [rw(i,1),pw(i,1)] = partialcorr(IFDw6_NGR(Vfd),Gf_cogs(Vfd,i),[DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
    [rw(i,2),pw(i,2)] = partialcorr(IFDw6_GSR(Vfd),Gf_cogs(Vfd,i),[DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);    
    [rw(i,3),pw(i,3)] = partialcorr(IFDw6_NGR(Vfd),Gf_cogs(Vfd,i),[DICV(Vfd),DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
    [rw(i,4),pw(i,4)] = partialcorr(IFDw6_GSR(Vfd),Gf_cogs(Vfd,i),[DICV(Vfd),DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
    [r_ICVxCog(i),p_ICVxCog(i)] = partialcorr(DICV(Vfd),Gf_cogs(Vfd,i),[DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
end

term_age = FixedEffect(DAge(Vfd), 'Age');
term_sex = FixedEffect(DSex(Vfd), 'Sex');
term_hd = FixedEffect(DHand(Vfd), 'Handed');
term_HM = FixedEffect(mFD(Vfd), 'mFD');
term_ICV = FixedEffect(DICV(Vfd), 'ICV');
term_FIQ = FixedEffect(Gf_cogs(Vfd,1), 'FIQ');
term_CF = FixedEffect(Gf_cogs(Vfd,3), 'CogFlex');

model_ICV = term_age+term_sex+term_hd+term_HM+term_ICV;
model_FIQ = term_age+term_sex+term_hd+term_HM+term_ICV+term_FIQ;
model_CF = term_age+term_sex+term_hd+term_HM+term_ICV+term_CF;


surfdata = {IFD_NGR_s8_data.cdata',IFD_GSR_s8_data.cdata'};

result_slm=cell(7,2);
for i = 1:2
    slm_ICV = SLM( model_ICV, DICV(Vfd), ...
        'surf', avsurf_hemi, ...
        'correction', {'rft', 'fdr'}, ...
        'cluster_threshold', 0.001, ...
        'mask', HemiMask,'two_tailed', false);
    slm_ICV.fit(surfdata{i}(Vfd,:));
    result_slm{1,i} = slm_ICV;


    for j=1:6
        term_cog = FixedEffect(Gf_cogs(Vfd,j), 'Cog');
        model_cog = term_age+term_sex+term_hd+term_HM+term_ICV+term_cog;

        slm_cog = SLM( model_cog, Gf_cogs(Vfd,j), ...
            'surf', avsurf_hemi, ...
            'correction', {'rft', 'fdr'}, ...
            'cluster_threshold', 0.001, ...
            'mask', HemiMask,'two_tailed', false);
        slm_cog.fit(surfdata{i}(Vfd,:));

        result_slm{1+j,i} = slm_cog;
    end

end


for i=1:7
    obj = plot_slm(result_slm{i,1}, avsurf_hemi, 'mask', HemiMask, 't_colorlimits',[-5,5]);
    pretty_plot(obj)
%     obj = plot_slm(result_slm{i,2}, avsurf_hemi, 'mask', HemiMask);
%     pretty_plot(obj)
end
%%
table_model = table(DICV(Vfd),DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd),IFDw6_NGR(Vfd));
glm1 = fitglm(table_model);%拟合方程
beta2 = glm1.Coefficients.Estimate(3);
beta3 = glm1.Coefficients.Estimate(4);
beta4 = glm1.Coefficients.Estimate(5);
beta5 = glm1.Coefficients.Estimate(6);
R = table2array(glm1.Residuals(:,2)); %残差
predict_fix = predict(glm1,table_model); 
f = predict_fix + R - beta2*table_model.Var2-beta3*table_model.Var3 - beta4*table_model.Var4 - beta5*table_model.Var5; 
x=DICV(Vfd);
save('ICVxIFD6_NGR.mat','x','f');

table_model = table(DICV(Vfd),DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd),Gf_cogs(Vfd,1));
glm1 = fitglm(table_model);%拟合方程
beta2 = glm1.Coefficients.Estimate(3);
beta3 = glm1.Coefficients.Estimate(4);
beta4 = glm1.Coefficients.Estimate(5);
beta5 = glm1.Coefficients.Estimate(6);
R = table2array(glm1.Residuals(:,2)); %残差
predict_fix = predict(glm1,table_model); 
f = predict_fix + R - beta2*table_model.Var2-beta3*table_model.Var3 - beta4*table_model.Var4 - beta5*table_model.Var5; 
x=DICV(Vfd);
save('ICVxFIQ.mat','x','f');

table_model = table(DICV(Vfd),DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd),Gf_cogs(Vfd,3));
glm1 = fitglm(table_model);%拟合方程
beta2 = glm1.Coefficients.Estimate(3);
beta3 = glm1.Coefficients.Estimate(4);
beta4 = glm1.Coefficients.Estimate(5);
beta5 = glm1.Coefficients.Estimate(6);
R = table2array(glm1.Residuals(:,2)); %残差
predict_fix = predict(glm1,table_model); 
f = predict_fix + R - beta2*table_model.Var2-beta3*table_model.Var3 - beta4*table_model.Var4 - beta5*table_model.Var5; 
x=DICV(Vfd);
save('ICVxCF.mat','x','f');

table_model = table(IFDw6_NGR(Vfd),DICV(Vfd),DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd),Gf_cogs(Vfd,1));
glm1 = fitglm(table_model);%拟合方程
beta2 = glm1.Coefficients.Estimate(3);
beta3 = glm1.Coefficients.Estimate(4);
beta4 = glm1.Coefficients.Estimate(5);
beta5 = glm1.Coefficients.Estimate(6);
beta6 = glm1.Coefficients.Estimate(7);
R = table2array(glm1.Residuals(:,2)); %残差
predict_fix = predict(glm1,table_model); 
f = predict_fix + R - beta2*table_model.Var2-beta3*table_model.Var3 - beta4*table_model.Var4 - beta5*table_model.Var5 - beta6*table_model.Var6; 
x=IFDw6_NGR(Vfd);
save('IFD6_NGRxFIQ.mat','x','f');

table_model = table(IFDw6_NGR(Vfd),DICV(Vfd),DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd),Gf_cogs(Vfd,3));
glm1 = fitglm(table_model);%拟合方程
beta2 = glm1.Coefficients.Estimate(3);
beta3 = glm1.Coefficients.Estimate(4);
beta4 = glm1.Coefficients.Estimate(5);
beta5 = glm1.Coefficients.Estimate(6);
beta6 = glm1.Coefficients.Estimate(7);
R = table2array(glm1.Residuals(:,2)); %残差
predict_fix = predict(glm1,table_model); 
f = predict_fix + R - beta2*table_model.Var2-beta3*table_model.Var3 - beta4*table_model.Var4 - beta5*table_model.Var5 - beta6*table_model.Var6; 
x=IFDw6_NGR(Vfd);
save('IFD6_NGRxCF.mat','x','f');

table_model = table(DAge(Vfd),DICV(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd),IFDw6_NGR(Vfd));
glm1 = fitglm(table_model);%拟合方程
beta2 = glm1.Coefficients.Estimate(3);
beta3 = glm1.Coefficients.Estimate(4);
beta4 = glm1.Coefficients.Estimate(5);
beta5 = glm1.Coefficients.Estimate(6);
R = table2array(glm1.Residuals(:,2)); %残差
predict_fix = predict(glm1,table_model); 
f = predict_fix + R - beta2*table_model.Var2-beta3*table_model.Var3 - beta4*table_model.Var4 - beta5*table_model.Var5; 
x=DAge(Vfd);
save('AgexIFD6_NGR.mat','x','f');

table_model = table(DHand(Vfd),DICV(Vfd),DSexN(Vfd),DAge(Vfd),mFD(Vfd),IFDw6_NGR(Vfd));
glm1 = fitglm(table_model);%拟合方程
beta2 = glm1.Coefficients.Estimate(3);
beta3 = glm1.Coefficients.Estimate(4);
beta4 = glm1.Coefficients.Estimate(5);
beta5 = glm1.Coefficients.Estimate(6);
R = table2array(glm1.Residuals(:,2)); %残差
predict_fix = predict(glm1,table_model); 
f = predict_fix + R - beta2*table_model.Var2-beta3*table_model.Var3 - beta4*table_model.Var4 - beta5*table_model.Var5; 
x=DHand(Vfd);
save('HandxIFD6_NGR.mat','x','f');

%%
table_model = table(DICV(Vfd),DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd),IFDw6_GSR(Vfd));
glm1 = fitglm(table_model);%拟合方程
beta2 = glm1.Coefficients.Estimate(3);
beta3 = glm1.Coefficients.Estimate(4);
beta4 = glm1.Coefficients.Estimate(5);
beta5 = glm1.Coefficients.Estimate(6);
R = table2array(glm1.Residuals(:,2)); %残差
predict_fix = predict(glm1,table_model); 
f = predict_fix + R - beta2*table_model.Var2-beta3*table_model.Var3 - beta4*table_model.Var4 - beta5*table_model.Var5; 
x=DICV(Vfd);
save('ICVxIFD6_GSR.mat','x','f');

table_model = table(IFDw6_GSR(Vfd),DICV(Vfd),DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd),Gf_cogs(Vfd,1));
glm1 = fitglm(table_model);%拟合方程
beta2 = glm1.Coefficients.Estimate(3);
beta3 = glm1.Coefficients.Estimate(4);
beta4 = glm1.Coefficients.Estimate(5);
beta5 = glm1.Coefficients.Estimate(6);
beta6 = glm1.Coefficients.Estimate(7);
R = table2array(glm1.Residuals(:,2)); %残差
predict_fix = predict(glm1,table_model); 
f = predict_fix + R - beta2*table_model.Var2-beta3*table_model.Var3 - beta4*table_model.Var4 - beta5*table_model.Var5 - beta6*table_model.Var6; 
x=IFDw6_NGR(Vfd);
save('IFD6_GSRxFIQ.mat','x','f');

%% global mediation

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


% vertex-wise mediation
vtsig = result_slm{2,1}.Q'<0.05;
sum(vtsig)
IFD_sig = mean(IFD_NGR_s8_data.cdata(vtsig,:))';
[r,p] = partialcorr(IFD_sig(Vfd),Gf_cogs(Vfd,1),[DICV(Vfd),DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
[r,p] = partialcorr(IFD_sig(Vfd),Gf_cogs(Vfd,1),[DICV(Vfd),DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd),IFDw6_NGR(Vfd)]);

X = zscore(DICV(Vfd));
Y = zscore(Gf_cogs(Vfd,1));
M = zscore(IFD_sig(Vfd));
covs = zscore([DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
[paths, stats] = mediation(X, Y, M, 'covs',covs, 'boot', 'verbose', 'bootsamples', 10000, 'doCIs');

vtsig = result_slm{4,1}.Q'<0.05;
IFD_sig = mean(IFD_NGR_s8_data.cdata(vtsig,:))';
[r,p] = partialcorr(IFD_sig(Vfd),Gf_cogs(Vfd,1),[DICV(Vfd),DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
[r,p] = partialcorr(IFD_sig(Vfd),Gf_cogs(Vfd,1),[DICV(Vfd),DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd),IFDw6_NGR(Vfd)]);

X = zscore(DICV(Vfd));
Y = zscore(Gf_cogs(Vfd,3));
M = zscore(IFD_sig(Vfd));
covs = zscore([DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
[paths, stats] = mediation(X, Y, M, 'covs',covs, 'boot', 'verbose', 'bootsamples', 10000, 'doCIs');
% mediation_path_diagram(stats)
% mediation_scatterplots(stats)


vtsig = result_slm{2,2}.Q'<0.05;
IFD_sig = mean(IFD_GSR_s8_data.cdata(vtsig,:))';
[r,p] = partialcorr(IFD_sig(Vfd),Gf_cogs(Vfd,1),[DICV(Vfd),DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);

X = zscore(DICV(Vfd));
Y = zscore(Gf_cogs(Vfd,1));
M = zscore(IFD_sig(Vfd));
covs = zscore([DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);
[paths, stats] = mediation(X, Y, M, 'covs',covs, 'boot', 'verbose', 'bootsamples', 10000, 'doCIs');

vtsig = result_slm{4,2}.Q'<0.05;
IFD_sig = mean(IFD_GSR_s8_data.cdata(vtsig,:))';
[r,p] = partialcorr(IFD_sig(Vfd),Gf_cogs(Vfd,1),[DICV(Vfd),DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);

X = zscore(DICV(Vfd));
Y = zscore(Gf_cogs(Vfd,3));
M = zscore(IFD_sig(Vfd));
covs = zscore([DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)]);

[paths, stats] = mediation(X, Y, M, 'covs',covs, 'boot', 'verbose', 'bootsamples', 10000, 'doCIs');
mediation_path_diagram(stats)
mediation_scatterplots(stats)


%%
GiftiTemp=gifti();
GiftiTemp.cdata=result_slm{2,1}.t';
save(GiftiTemp,'/Volumes/eSSD/FinalCodes/ResultsMaps/IFDxCog/IFDxFIQ_Tmap_NGR.shape.gii','Base64Binary');
GiftiTemp.cdata=result_slm{2,1}.Q';
save(GiftiTemp,'/Volumes/eSSD/FinalCodes/ResultsMaps/IFDxCog/IFDxFIQ_Pmap_FDR_NGR.shape.gii','Base64Binary');
GiftiTemp.cdata=result_slm{4,1}.t';
save(GiftiTemp,'/Volumes/eSSD/FinalCodes/ResultsMaps/IFDxCog/IFDxDSST_Tmap_NGR.shape.gii','Base64Binary');
GiftiTemp.cdata=result_slm{4,1}.Q';
save(GiftiTemp,'/Volumes/eSSD/FinalCodes/ResultsMaps/IFDxCog/IFDxDSST_Pmap_FDR_NGR.shape.gii','Base64Binary');

GiftiTemp.cdata=result_slm{2,2}.t';
save(GiftiTemp,'/Volumes/eSSD/FinalCodes/ResultsMaps/IFDxCog/IFDxFIQ_Tmap_GSR.shape.gii','Base64Binary');
GiftiTemp.cdata=result_slm{2,2}.Q';
save(GiftiTemp,'/Volumes/eSSD/FinalCodes/ResultsMaps/IFDxCog/IFDxFIQ_Pmap_FDR_GSR.shape.gii','Base64Binary');
GiftiTemp.cdata=result_slm{4,2}.t';
save(GiftiTemp,'/Volumes/eSSD/FinalCodes/ResultsMaps/IFDxCog/IFDxDSST_Tmap_GSR.shape.gii','Base64Binary');
GiftiTemp.cdata=result_slm{4,2}.Q';
save(GiftiTemp,'/Volumes/eSSD/FinalCodes/ResultsMaps/IFDxCog/IFDxDSST_Pmap_FDR_GSR.shape.gii','Base64Binary');

GiftiTemp=gifti();
GiftiTemp.cdata=result_slm{2,1}.t';
save(GiftiTemp,'/Volumes/eSSD/FinalCodes/ResultsMaps/IFDxCog/IFDxFIQ_Uri387_Tmap_NGR.shape.gii','Base64Binary');
GiftiTemp.cdata=result_slm{2,1}.Q';
save(GiftiTemp,'/Volumes/eSSD/FinalCodes/ResultsMaps/IFDxCog/IFDxFIQ_Uri387_Pmap_FDR_NGR.shape.gii','Base64Binary');
GiftiTemp.cdata=result_slm{4,1}.t';
save(GiftiTemp,'/Volumes/eSSD/FinalCodes/ResultsMaps/IFDxCog/IFDxDSST_Uri387_Tmap_NGR.shape.gii','Base64Binary');
GiftiTemp.cdata=result_slm{4,1}.Q';
save(GiftiTemp,'/Volumes/eSSD/FinalCodes/ResultsMaps/IFDxCog/IFDxDSST_Uri387_Pmap_FDR_NGR.shape.gii','Base64Binary');
