function obj = plot_slm(slm, surfaces, options)
% Tutorial function for plotting SLM data.
%
%   obj = TUTORIAL_PLOT_SLM(slm, surfaces, options) plots the data in an slm
%   to the cortical surfaces included in cell array surfaces. The following
%   name-value pairs are allowed:
%
%   mask
%       A logical vector with `false` for all vertices to exclude.
%   plot_t 
%       Plot the t-values, defaults to true.
%   plot_clus
%       Plot the clusterwise p-values, defaults to false. 
%   plot_peak
%       Plot the peak p-values, defaults to false.
%   plot_fdr
%       Plot the vertexwise p-values, defaults to false.
%   alpha
%       The upper limit of p-values to plot, defaults to 0.05.
%   t_colorlimits
%       The lower/upper limit of t-values to plot, defaults to [-6, 4].

arguments
    slm
    surfaces
    options.mask (:,1) = nan
    options.alpha (1,1) double = 0.05
    options.t_colorlimits (1,2) double = [-4, 4]
end

to_plot = [];
labels = {};
colormaps = {};
colorlimits = [];

to_plot = [to_plot, slm.t(:)];
labels = [labels, {'t-values'}];
colormaps = [colormaps; {[parula; .7 .7 .7]}];
colorlimits = [colorlimits; options.t_colorlimits];

% tempclus= ones(1,length(slm.P.clusid{1,1}));
% for i = 1:max(slm.P.clus{1,1}.clusid)
%     tempclus(slm.P.clusid{1,1}==i) = slm.P.clus{1,1}.P(i);
% end
% to_plot = [to_plot, tempclus(:)];
% labels = [labels, {{'Cluster', 'p-values (RFT)'}}];
% colormaps = [colormaps; {[flipud(autumn); .7 .7 .7]}];
% colorlimits = [colorlimits; 0, options.alpha];

to_plot = [to_plot, slm.Q(:)];
labels = [labels, {{'p-values (FDR)'}}];
colormaps = [colormaps; {[flipud(autumn); .7 .7 .7]}];
colorlimits = [colorlimits; 0, 0.05];

if ~isnan(options.mask)
    to_plot(~options.mask, :) = inf;
end

obj = plot_hemispheres(...
    to_plot,  ...
    surfaces, ...
    'labeltext', labels ...
    );

obj.colormaps(colormaps);
obj.colorlimits(colorlimits);
end