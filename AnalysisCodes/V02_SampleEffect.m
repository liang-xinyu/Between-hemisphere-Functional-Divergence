% load data (Distance in diverse dimensions)
Dists=load('/Volumes/eSSD/FinalCodes/CalcData/ALL777_IFD.mat');
IFDw6_NGR=Dists.IFdistAll_AVG{1}(:,6); %NGR
IFDw6_GSR=Dists.IFdistAll_AVG{2}(:,6); %GSR
conID=Dists.SubID;
clear Dists

%% varis
SubSize=length(conID);
HeadMotion = load('/Volumes/eSSD/Project-HCP-HemiSurfAsym/GradientAna/HeadMotion_FD_777subs.mat');
[~,~,ihm]=intersect(conID,HeadMotion.conID,'stable');
mFD = (HeadMotion.mFD1(ihm)+HeadMotion.mFD2(ihm))/2;

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


Xvar = IFDw6_NGR(Vfd);
Yvar = Gf_cogs(Vfd,1);
Covar = [DAge(Vfd),DSexN(Vfd),DHand(Vfd),mFD(Vfd)];

R_resample = zeros(11,1000);
P_resample = zeros(11,1000);
for i =1:1000
    s25 = randsample([1:755],25);
    s33 = randsample([1:755],33);
    s50 = randsample([1:755],50);
    s70 = randsample([1:755],70);
    s100 = randsample([1:755],100);
    s135 = randsample([1:755],135);
    s200 = randsample([1:755],200);
    s265 = randsample([1:755],265);
    s375 = randsample([1:755],375);
    s525 = randsample([1:755],525);
    s725 = randsample([1:755],725);
    [R_resample(1,i), P_resample(1,i)]= partialcorr(Xvar(s25),Yvar(s25),Covar(s25,:));
    [R_resample(2,i), P_resample(2,i)]= partialcorr(Xvar(s33),Yvar(s33),Covar(s33,:));
    [R_resample(3,i), P_resample(3,i)]= partialcorr(Xvar(s50),Yvar(s50),Covar(s50,:));
    [R_resample(4,i), P_resample(4,i)]= partialcorr(Xvar(s70),Yvar(s70),Covar(s70,:));
    [R_resample(5,i), P_resample(5,i)]= partialcorr(Xvar(s100),Yvar(s100),Covar(s100,:));
    [R_resample(6,i), P_resample(6,i)]= partialcorr(Xvar(s135),Yvar(s135),Covar(s135,:));
    [R_resample(7,i), P_resample(7,i)]= partialcorr(Xvar(s200),Yvar(s200),Covar(s200,:));
    [R_resample(8,i), P_resample(8,i)]= partialcorr(Xvar(s265),Yvar(s265),Covar(s265,:));
    [R_resample(9,i), P_resample(9,i)]= partialcorr(Xvar(s375),Yvar(s375),Covar(s375,:));
    [R_resample(10,i), P_resample(10,i)]= partialcorr(Xvar(s525),Yvar(s525),Covar(s525,:));
    [R_resample(11,i), P_resample(11,i)]= partialcorr(Xvar(s725),Yvar(s725),Covar(s725,:));
end



R_resample = zeros(15,1000);
P_resample = zeros(15,1000);
for i =1:1000
    s50 = randsample([1:755],50);
    s100 = randsample([1:755],100);
    s150 = randsample([1:755],150);
    s200 = randsample([1:755],200);
    s250 = randsample([1:755],250);
    s300 = randsample([1:755],300);
    s350 = randsample([1:755],350);
    s400 = randsample([1:755],400);
    s450 = randsample([1:755],450);
    s500 = randsample([1:755],500);
    s550 = randsample([1:755],550);
    s600 = randsample([1:755],600);
    s650 = randsample([1:755],650);
    s700 = randsample([1:755],700);
    s750 = randsample([1:755],750);
    [R_resample(1,i), P_resample(1,i)]= partialcorr(Xvar(s50),Yvar(s50),Covar(s50,:));
    [R_resample(2,i), P_resample(2,i)]= partialcorr(Xvar(s100),Yvar(s100),Covar(s100,:));
    [R_resample(3,i), P_resample(3,i)]= partialcorr(Xvar(s150),Yvar(s150),Covar(s150,:));
    [R_resample(4,i), P_resample(4,i)]= partialcorr(Xvar(s200),Yvar(s200),Covar(s200,:));
    [R_resample(5,i), P_resample(5,i)]= partialcorr(Xvar(s250),Yvar(s250),Covar(s250,:));
    [R_resample(6,i), P_resample(6,i)]= partialcorr(Xvar(s300),Yvar(s300),Covar(s300,:));
    [R_resample(7,i), P_resample(7,i)]= partialcorr(Xvar(s350),Yvar(s350),Covar(s350,:));
    [R_resample(8,i), P_resample(8,i)]= partialcorr(Xvar(s400),Yvar(s400),Covar(s400,:));
    [R_resample(9,i), P_resample(9,i)]= partialcorr(Xvar(s450),Yvar(s450),Covar(s450,:));
    [R_resample(10,i), P_resample(10,i)]= partialcorr(Xvar(s500),Yvar(s500),Covar(s500,:));
    [R_resample(11,i), P_resample(11,i)]= partialcorr(Xvar(s550),Yvar(s550),Covar(s550,:));
    [R_resample(12,i), P_resample(12,i)]= partialcorr(Xvar(s600),Yvar(s600),Covar(s600,:));
    [R_resample(13,i), P_resample(13,i)]= partialcorr(Xvar(s650),Yvar(s650),Covar(s650,:));
    [R_resample(14,i), P_resample(14,i)]= partialcorr(Xvar(s700),Yvar(s700),Covar(s700,:));
    [R_resample(15,i), P_resample(15,i)]= partialcorr(Xvar(s750),Yvar(s750),Covar(s750,:));
end

R_inv=R_resample';
P_inv=P_resample';