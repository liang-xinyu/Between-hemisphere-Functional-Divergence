% Generate FGD

R1_NGR_path=g_ls('/Volumes/eSSD/FinalCodes/CalcData/Realigned_ALL777_Com/R1NGR/*.mat');
R2_NGR_path=g_ls('/Volumes/eSSD/FinalCodes/CalcData/Realigned_ALL777_Com/R2NGR/*.mat');
R1_GSR_path=g_ls('/Volumes/eSSD/FinalCodes/CalcData/Realigned_ALL777_Com/R1GSR/*.mat');
R2_GSR_path=g_ls('/Volumes/eSSD/FinalCodes/CalcData/Realigned_ALL777_Com/R2GSR/*.mat');
IDfile=load('/Volumes/eSSD/FinalCodes/Inform/Subinfo_rhand_nocode.mat');
conID=IDfile.Sub_ID;
savepath='/Volumes/eSSD/FinalCodes/CalcData/';
IFD_calc_save(savepath,'ALL777_IFD',conID,R1_NGR_path,R2_NGR_path,R1_GSR_path,R2_GSR_path)

clear
R1_NGR_path=g_ls('/Volumes/eSSD/FinalCodes/CalcData/Realigned_387_Com/R1NGR/*.mat');
R2_NGR_path=g_ls('/Volumes/eSSD/FinalCodes/CalcData/Realigned_387_Com/R2NGR/*.mat');
R1_GSR_path=g_ls('/Volumes/eSSD/FinalCodes/CalcData/Realigned_387_Com/R1GSR/*.mat');
R2_GSR_path=g_ls('/Volumes/eSSD/FinalCodes/CalcData/Realigned_387_Com/R2GSR/*.mat');
IDfile=load('/Volumes/eSSD/FinalCodes/Inform/Sub_UnrelatedID.mat');
conID=IDfile.Unrelated_A;
savepath='/Volumes/eSSD/FinalCodes/CalcData/';
IFD_calc_save(savepath,'Uri387_IFD',conID,R1_NGR_path,R2_NGR_path,R1_GSR_path,R2_GSR_path)

%
clear
R1_NGR_path=g_ls('/Volumes/eSSD/FinalCodes/CalcData/Residual366_to387Com/R1NGR/*.mat');
R2_NGR_path=g_ls('/Volumes/eSSD/FinalCodes/CalcData/Residual366_to387Com/R2NGR/*.mat');
R1_GSR_path=g_ls('/Volumes/eSSD/FinalCodes/CalcData/Residual366_to387Com/R1GSR/*.mat');
R2_GSR_path=g_ls('/Volumes/eSSD/FinalCodes/CalcData/Residual366_to387Com/R2GSR/*.mat');
IDfile=load('/Volumes/eSSD/FinalCodes/Inform/Sub_UnrelatedID.mat');
Unrid=[IDfile.Unrelated_B;IDfile.Unrelated_C];
conID=sort(Unrid);
savepath='/Volumes/eSSD/FinalCodes/CalcData/';
IFD_calc_save(savepath,'Res366_IFD',conID,R1_NGR_path,R2_NGR_path,R1_GSR_path,R2_GSR_path)

%% Generate FGD based on R2L
clear
R1_NGR_path=g_ls('/Volumes/eSSD/FinalCodes/RealignedRtoL_ALL777/R1NGR/*.mat');
R2_NGR_path=g_ls('/Volumes/eSSD/FinalCodes/RealignedRtoL_ALL777/R2NGR/*.mat');
R1_GSR_path=g_ls('/Volumes/eSSD/FinalCodes/RealignedRtoL_ALL777/R1GSR/*.mat');
R2_GSR_path=g_ls('/Volumes/eSSD/FinalCodes/RealignedRtoL_ALL777/R2GSR/*.mat');
IDfile=load('/Volumes/eSSD/FinalCodes/Inform/Subinfo_rhand_nocode.mat');
conID=IDfile.Sub_ID;
savepath='/Volumes/eSSD/FinalCodes/CalcData/';
IFD_calc_save(savepath,'ALL777_R2L_IFD',conID,R1_NGR_path,R2_NGR_path,R1_GSR_path,R2_GSR_path)
