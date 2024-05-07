% Please change this path according to your saved path 
Codes_path = '/Volumes/eSSD/Between-hemisphere-Functional-Divergence/';

R1=load([Codes_path,'Data/TemplateGradients/TemplateMat/GroupAVG_777_R1NGR_gradient.mat']);
R2=load([Codes_path,'Data/TemplateGradients/TemplateMat/GroupAVG_777_R2NGR_gradient.mat']);

explain_R1=R1.g_hemi.lambda{1}(1:10)/sum(R1.g_hemi.lambda{1}(1:10));
explain_R2=R2.g_hemi.lambda{1}(1:10)/sum(R2.g_hemi.lambda{1}(1:10));

explain_R1_all=R1.g_hemi.lambda{1}/sum(R1.g_hemi.lambda{1});
explain_R2_all=R2.g_hemi.lambda{1}/sum(R2.g_hemi.lambda{1});

% This file need to be download separately 
Dists=load([Codes_path,'Data/ALL777_IFD.mat'],'ICC_IFdist_whole');
ICC=Dists.ICC_IFdist_whole(:,1);

explain10=(explain_R1+explain_R2)/2;
explainAlL=(explain_R1_all+explain_R2_all)/2;
