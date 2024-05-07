function [smoothedpath] = SurfSmooth_wb(Vari4s,kernel,HemiMask,outputname)
%% smoothing using wb_command 
% Vari4s, the cortical map for smoothing
% Kernel, the size of smoothing kernel (fwhm)
% HemiMask, medial wall mask
% outputname, file name to save as 

% filepath for saving files
workfolder='/Volumes/eSSD/FinalCodes/Folder4Smooth/';

metric=zeros(size(Vari4s,1),10242);
metric(:,HemiMask)=Vari4s;

% load a gifti file as template for saving
template=gifti('/Volumes/eSSD/FinalCodes/ResultsMaps/IFD777_NGR_AVG_left.shape.gii');
template.cdata=metric';
metricpath=[workfolder,filesep,outputname,'_wait4smooth.shape.gii'];
save(template,metricpath,'Base64Binary');

% smoothing based on the left hemisphere
surfLpath='/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/S1200.L.midthickness_MSMAll.10k_fs_LR.surf.gii';
% set up wb_command
wbc='/Applications/workbench/bin_macosx64/wb_command';
% filepath of output
smoothedpath=[workfolder,filesep,outputname,'_smoothed_s',num2str(kernel),'.shape.gii'];

smoothcmd = [wbc,' -metric-smoothing ',surfLpath,' ',metricpath,' ',num2str(kernel),' ',smoothedpath,...
    ' -fwhm -fix-zeros'];
system(smoothcmd)

end

