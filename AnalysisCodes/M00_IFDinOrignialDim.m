mwL=gifti('D:\Project-HCP\HCP_Atlas\Atlas-fs_LR_10k\fs_LR_10k\Human.MedialWall_Conte69.L.10k_fs_LR.label.gii');
mwR=gifti('D:\Project-HCP\HCP_Atlas\Atlas-fs_LR_10k\fs_LR_10k\Human.MedialWall_Conte69.R.10k_fs_LR.label.gii');
HemiMask=logical((1-mwL.cdata') .* (1-mwR.cdata')); 

surfL=gifti('D:\Project-HCP\HCP_Atlas\Atlas-fs_LR_10k\fs_LR_10k\S1200.L.inflated_MSMAll.10k_fs_LR.surf.gii');
avsurf_hemi.coord=surfL.vertices';
avsurf_hemi.tri=surfL.faces;

NGR_L=(R1NGR_L+R2NGR_L)/2;
NGR_R=(R1NGR_R+R2NGR_R)/2;
clear R1NGR_L R2NGR_L R1NGR_R R2NGR_R
NGR_L=NGR_L(HemiMask,HemiMask);
NGR_R=NGR_R(HemiMask,HemiMask);
NGR_hemi=(NGR_L+NGR_R)/2;

NGR_hemi_cs = 1-squareform(pdist(NGR_hemi','cosine'));
NGR_hemi_na = 1-acos(NGR_hemi_cs)/pi;

NGR_L_cs = 1-squareform(pdist(NGR_L','cosine'));
NGR_L_na = 1-acos(NGR_L_cs)/pi;

NGR_R_cs = 1-squareform(pdist(NGR_R','cosine'));
NGR_R_na = 1-acos(NGR_R_cs)/pi;

[realigned_L, ~] = procrustes_analysis({NGR_L_na},10,NGR_hemi_na);
[realigned_R, ~] = procrustes_analysis({NGR_R_na},10,NGR_hemi_na);

dist_highdim=zeros(9354,1);
dist_highdim_orin=zeros(9354,1);
dist_FC=zeros(9354,1);
for i=1:9354
    dist_highdim(i,1)=dist(realigned_L(i,:),realigned_R(i,:)');
    dist_highdim_orin(i,1)=dist(NGR_L_na(i,:),NGR_R_na(i,:)');
end

for i=1:9354
    dist_FC(i,1)=dist(NGR_L(i,:),NGR_R(i,:)');
end


hemigra=zeros(1,10242);
hemigra(HemiMask)=dist_highdim;
figure,SurfStatViewData_lxy(hemigra,avsurf_hemi,[0,6],'','Hori');

hemigra=zeros(1,10242);
hemigra(HemiMask)=dist_highdim_orin;
figure,SurfStatViewData_lxy(hemigra,avsurf_hemi,[0,6],'','Hori');

hemigra=zeros(1,10242);
hemigra(HemiMask)=dist_FC;
figure,SurfStatViewData_lxy(hemigra,avsurf_hemi,[0,6],'','Hori');

[rld_L, ~] = procrustes_analysis(g_hemi_L.gradients,10,g_hemi.gradients{1});
[rld_R, ~] = procrustes_analysis(g_hemi_R.gradients,10,g_hemi.gradients{1});

dist_gradient=zeros(9354,1);
for i=1:9354
    dist_gradient(i,1)=dist(rld_L(i,:),rld_R(i,:)');
end
hemigra=zeros(1,10242);
hemigra(HemiMask)=dist_gradient;
figure,SurfStatViewData_lxy(hemigra,avsurf_hemi,[0,0.15],'','Hori');

%%
dist_FC_100206=zeros(9354,1);
hcon_L_masked=hcon_L(HemiMask,HemiMask);
hcon_R_masked=hcon_R(HemiMask,HemiMask);
for i=1:9354
    dist_FC_100206(i,1)=dist(hcon_L_masked(i,:),hcon_R_masked(i,:)');
end
hemigra=zeros(1,10242);
hemigra(HemiMask)=dist_FC_100206;
figure,SurfStatViewData_lxy(hemigra,avsurf_hemi,[5,20],'','Hori');

[rld_L_100206, ~] = procrustes_analysis(g_hemiL.gradients,10,R1NGR_Grads);
[rld_R_100206, ~] = procrustes_analysis(g_hemiR.gradients,10,R1NGR_Grads);

dist_gradient_100206=zeros(9354,1);
for i=1:9354
    dist_gradient_100206(i,1)=dist(rld_L_100206(i,:),rld_R_100206(i,:)');
end
hemigra=zeros(1,10242);
hemigra(HemiMask)=dist_gradient_100206;
figure,SurfStatViewData_lxy(hemigra,avsurf_hemi,[0,0.15],'','Hori');

[prctile_rank,r_original,r_rand]=spintest_10k(dist_gradient_100206,dist_FC_100206,HemiMask,5000,'spin');

intraFCS_L_100206=sum(hcon_L);
intraFCS_R_100206=sum(hcon_R);
intraFCSpos_L_100206=intraFCS_L_100206.*(intraFCS_L_100206>0);
intraFCSpos_R_100206=intraFCS_R_100206.*(intraFCS_R_100206>0);

LI_iFCS_100206=abs(intraFCSpos_L_100206-intraFCSpos_R_100206)./(intraFCSpos_L_100206+intraFCSpos_R_100206);
LI_iFCS_100206(isnan(LI_iFCS_100206))=0;

dist_iFCS_100206=abs(intraFCS_L_100206-intraFCS_R_100206);
figure,SurfStatViewData_lxy(dist_iFCS_100206,avsurf_hemi,[0,1300],'','Hori');
figure,SurfStatViewData_lxy(LI_iFCS_100206,avsurf_hemi,[0,1],'','Hori');

SurfStatPlot(dist_gradient_100206,dist_iFCS_100206(HemiMask)')
[prctile_rank,r_original,r_rand]=spintest_10k(dist_gradient_100206,dist_iFCS_100206(HemiMask)',HemiMask,5000,'spin');
