function [prctile_rank,r_original,r_rand]=spintest_10k(X,Y,HemiMask,npt,method)

if strcmp(method,'spin') % spin permutation
    % load sphere
    spL10k=gifti('/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/Sphere.10k.L.surf.gii');
    sphereL.coord=spL10k.vertices';
    sphereL.tri=spL10k.faces;
    
    for k=1:size(X,2)
        % Template functional gradient
        Xfull=zeros(1,10242);
        Yfull=zeros(1,10242);
        Xfull(HemiMask)=X(:,k);
        Yfull(HemiMask)=Y;
        
        % Let's create some rotations
        n_permutations = npt;
        tic
        y_rand = spin_permutations({Yfull'},{sphereL},n_permutations,'random_state',0);
        toc
        
        % Merge the rotated data into single vectors
        Y_rotated = squeeze(y_rand{1}(:,1,:));
        
        % Find correlation between FC-G1 with thickness and T1w/T2w
        flag_full=(Xfull~=0)&(Yfull~=0);
        [r_original, ~] = corr(Xfull(flag_full)',Yfull(flag_full)','rows','pairwise','type','spearman');
        r_rand=zeros(1,n_permutations);
        for j=1:n_permutations
             flag_perm=(Xfull~=0)&(Y_rotated(:,j)~=0)';
             r_rand(1,j) = corr(Xfull(flag_perm)',Y_rotated(flag_perm',j),'rows','pairwise','type','spearman');
        end
        prctile_rank(k) = mean(r_original > r_rand);
        %significant = prctile_rank < 0.025 || prctile_rank >= 0.975;

    end
elseif strcmp(method,'moran') % Moran Spectral Randomization (MSR)
    % load surface
    surfL=gifti('/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/S1200.L.midthickness_MSMAll.10k_fs_LR.surf.gii');
    avsurf_hemi.coord=surfL.vertices';
    avsurf_hemi.tri=surfL.faces;
            
    n_ring = 1;
    tic
    MEM = compute_mem(avsurf_hemi,'n_ring',n_ring,'mask',~HemiMask);
    toc
    for k=1:size(X,2)
        n_rand = npt;
        tic
        y_rand = moran_randomization(Y',MEM,n_rand,'procedure','singleton','joint',true,'random_state',0);
        Y_spinned = squeeze(y_rand(:,1,:));
        toc
        [r_original(k,1), ~] = corr(X(:,k),Y','rows','pairwise','type','pearson');
        r_rand(k,:) = corr(X(:,k),Y_spinned,'rows','pairwise','type','pearson');
        prctile_rank(k) = mean(r_original(k,1) > r_rand(k,:));
    end

elseif strcmp(method,'surrogate')
    surfL=gifti('/Volumes/eSSD/Accessory/HCP_Atlas/Atlas-fs_LR_10k/fs_LR_10k/S1200.L.midthickness_MSMAll.10k_fs_LR.surf.gii');
    avsurf_hemi.coord=surfL.vertices';
    avsurf_hemi.tri=surfL.faces;

    G = surface_to_graph(avsurf_hemi, 'geodesic', ~HemiMask, true);
    geodesic_distance = distances(G);

    random_initialization = 0;
    n_surrogate_datasets = 1000;
    
    % Note: num_samples must be greater than num_neighbors
    num_samples = inf;
    num_neighbors = 100;
    
    obj_subsample = variogram(geodesic_distance, 'ns', num_samples, ...
        'knn', num_neighbors, 'random_state', random_initialization, 'num_workers',4);
    tic
    surrogates_subsample = obj_subsample.fit(Y, n_surrogate_datasets);
    toc

    r_real = corr(Y, X);
    r_surrogate = corr(surrogates_subsample, X);
    prctile_rank = mean(r_real > r_surrogate);
end


% Compute percentile rank.
%  significant_pdist = prctile_rank_pdist < 0.025 || prctile_rank_pdist >= 0.975;
%p_spin = 1-prctile_rank;
%  significant_pdist = prctile_rank_pdist < 0.025 || prctile_rank_pdist >= 0.975;

