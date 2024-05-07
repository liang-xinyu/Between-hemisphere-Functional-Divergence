% This is a demo code to show how we measured the between-hemisphere functional distance
% The original files are too large, so they are not uploaded

% Please change this path according to your saved path 
Codes_path = '/Volumes/eSSD/Between-hemisphere-Functional-Divergence/';

% load files includes the preprocessed and realigned gradients with
% different conditions

R1_NGR_path=g_ls([Codes_path,'CalcData/Realigned_ALL777_Com/R1NGR/*.mat']);
R2_NGR_path=g_ls([Codes_path,'CalcData/Realigned_ALL777_Com/R2NGR/*.mat']);
R1_GSR_path=g_ls([Codes_path,'CalcData/Realigned_ALL777_Com/R1GSR/*.mat']);
R2_GSR_path=g_ls([Codes_path,'CalcData/Realigned_ALL777_Com/R2GSR/*.mat']);

% Selected subjects in HCP
IDfile=load([Codes_path,'Data/Subinfo_rhand_nocode.mat']);
conID=IDfile.Sub_ID;
savepath=[Codes_path,'Data/CalcData/'];

% Call function and generated between-hemisphere functional distance (denote as IFD[interhemispheric functional distance] in codes)
IFD_calc_save(savepath,'ALL777_IFD',conID,R1_NGR_path,R2_NGR_path,R1_GSR_path,R2_GSR_path)


