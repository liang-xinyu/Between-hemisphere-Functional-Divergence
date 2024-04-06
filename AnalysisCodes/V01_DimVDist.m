IFD=load('/Volumes/eSSD/FinalCodes/CalcData/ALL777_IFD.mat');
IFdist6_NGR=IFD.IFdist_AVG{6,1}; %NGR
IFdist6_ID = IFD.SubID;
%IFdist6_GSR=IFD.IFdist_AVG{6,2}; %GSR
clear IFD

R1_NGR_path=g_ls('/Volumes/eSSD/FinalCodes/CalcData/Realigned_ALL777_Com/R1NGR/*.mat');
R2_NGR_path=g_ls('/Volumes/eSSD/FinalCodes/CalcData/Realigned_ALL777_Com/R2NGR/*.mat');

IDfile=load('/Volumes/eSSD/FinalCodes/Inform/Subinfo_rhand_nocode.mat');
conID=IDfile.Sub_ID;

SubSize=length(conID);
gmats=zeros(SubSize,9354);
IFdiff=repmat({gmats,gmats},6,1);
for i=1:SubSize
    tmpg1=load(R1_NGR_path{i});
    tmpg2=load(R2_NGR_path{i});   
    for j =1:6
        gradL1=tmpg1.realigned_L(:,j);
        gradR1=tmpg1.realigned_R(:,j);
        IFdiff{j,1}(i,:)=abs(gradL1 - gradR1); % R1 NGR
        
        gradL2=tmpg2.realigned_L(:,j);
        gradR2=tmpg2.realigned_R(:,j);
        IFdiff{j,2}(i,:)=abs(gradL2 - gradR2); % R2 NGR
    end 
end

diffsimi = zeros(777,6);
for subj = 1:777
    for w = 1:6
        meandiff= (IFdiff{w,1}(subj,:) + IFdiff{w,2}(subj,:))/2;
        diffsimi(subj,w) = corr(IFdist6_NGR(subj,:)',meandiff');
    end
end

diffsimiind = zeros(6,1);
diffsimiind_p = zeros(6,1);

for w =1:6
    meandiff= (IFdiff{w,1} + IFdiff{w,2})/2;
    [diffsimiind(w,:),diffsimiind_p(w,:)] = corr(mean(IFdist6_NGR,2),mean(meandiff,2));
end

save('diffsime_dimensions.mat','diffsimi','diffsimiind');
