function [ a, cb ] = SurfStatViewData_lxy( data, surf, caxis_range, figtitle, BrainType)

%Basic viewer for surface data.
% 
% Usage: [ a, cb ] = SurfStatViewData( data, surf [,title [,background]] );
% 
% data        = 1 x v vector of data, v=#vertices
% surf.coord  = 3 x v matrix of coordinates.
% surf.tri    = t x 3 matrix of triangle indices, 1-based, t=#triangles.
% title       = any string, data name by default.
% background  = background colour, any matlab ColorSpec, such as 
%   'white' (default), 'black'=='k', 'r'==[1 0 0], [1 0.4 0.6] (pink) etc.
%   Letter and line colours are inverted if background is dark (mean<0.5).
%
% a  = vector of handles to the axes, left to right, top to bottom. 
% cb = handle to the colorbar.

% BrainType='Vert';
% BrainType='Hori';


% find cut between hemispheres, assuming they are concatenated
t=size(surf.tri,1);
v=size(surf.coord,2);
tmax=max(surf.tri,[],2);
tmin=min(surf.tri,[],2);

% to save time, check that the cut is half waystatspath
if min(tmin(t/2+1:t))-max(tmax(1:t/2))==1
    cut=t/2;
    cuv=v/2;
else % check all cuts
    for i=1:t-1
        tmax(i+1)=max(tmax(i+1),tmax(i));
        tmin(t-i)=min(tmin(t-i),tmin(t-i+1));
    end
    cut=min([find((tmin(2:t)-tmax(1:t-1))==1) t]);
    cuv=tmax(cut);
end

tl=1:cut;
tr=(cut+1):t;
vl=1:cuv;
vr=(cuv+1):v;

clf;

%% colormap is hot
% cm=[ones(1,3)*0.7;
%     ones(511,1)    (0:510)'/511  zeros(511,1)];
% cm(1,:)=ones(1,3)*0.7;

%% colormap is spectral
% cm=spectral(512);
% cm(end-26:end,:)=[];
% cm(1:75,:)=[];

%% colormap is parula
% cm=parula;
% 
% % cmorig=parula;
% % n = length(cmorig);
% % X0 = linspace (1, n, 1024);
% % cm = interp1(1:n,cmorig,X0);

%% colormap is jet
cmorig=jet;
n = length(cmorig);
X0 = linspace (1, n, 1024);
cm = interp1(1:n,cmorig,X0);

%% set zero to gray

clim = [min(caxis_range),max(caxis_range)];
if clim(1) ==0
%     cm=cm(round(size(cm,1)/2):end,:);
    Pos=[1];
    cm(Pos, :)=repmat([0.75, 0.75, 0.75],1,1);
elseif clim(2) ==0
%     cm=cm(1:round(size(cm,1)/2),:);
    Pos=[size(cm,1)-1:size(cm,1)];
    cm(Pos, :)=repmat([0.75, 0.75, 0.75],2,1);
   
elseif clim(1)*clim(2) > 0
    Pos=0;  
elseif clim(1)+clim(2)==0
    Pos=floor(size(cm,1)/2)+1;
    cm(Pos, :)=[0.75, 0.75, 0.75];
else   
    Pos=floor(size(cm,1)*abs(clim(1))/(abs(clim(1))+abs(clim(2))))+1;
    cm(Pos, :)=[0.75, 0.75, 0.75];

end

colormap(cm);
% 
%         cmp=colormap(gca);
% 
%         colormap(gca, cmp);
%         
% clim=[min(data),max(data)];
% if clim(1)==clim(2)
%     clim=clim(1)+[-1 0];
% end

if cut<t
    
    % backgroud size
    tmp=get(gcf,'Position');
    set(gcf,'Position', [50, tmp(2) ,tmp(4)*2.4, tmp(4)*1.8]);
    
    % surface position and parameters
    % left dorsal
    a(1)=axes('position',[0.05 0.24 0.4 1]);
    trisurf(surf.tri(tl,:),surf.coord(1,vl),surf.coord(2,vl),surf.coord(3,vl),...
        double(data(vl)),'EdgeColor','none');
    view(-90,0);
    daspect([1 1 1]); axis tight; camlight; axis vis3d off; zoom(0.95);
    lighting phong; material([0.5 0.5 0.2 0.2]); shading flat;
    caxis(caxis_range)
    
    % left ventral
    a(3)=axes('position',[0.05 -0.24 0.4 1]);
    trisurf(surf.tri(tl,:),surf.coord(1,vl),surf.coord(2,vl),surf.coord(3,vl),...
        double(data(vl)),'EdgeColor','none');
    view(90,0);
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;zoom(0.95);
    lighting phong; material([0.5 0.5 0.2 0.2]); shading flat;
    caxis(caxis_range)
    
    % right dorsal
    a(2)=axes('position',[0.48 0.24 0.4 1]);
    trisurf(surf.tri(tr,:)-cuv,surf.coord(1,vr),surf.coord(2,vr),surf.coord(3,vr),...
        double(data(vr)),'EdgeColor','none');
    view(90,0);
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;zoom(0.95);
    lighting phong; material([0.5 0.5 0.2 0.2]); shading flat;
    caxis(caxis_range)
    
    % right ventral
    a(4)=axes('position',[0.48 -0.24 0.4 1]);
    trisurf(surf.tri(tr,:)-cuv,surf.coord(1,vr),surf.coord(2,vr),surf.coord(3,vr),...
        double(data(vr)),'EdgeColor','none');
    view(-90,0);
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;zoom(0.95);
    lighting phong; material([0.5 0.5 0.2 0.2]); shading flat;
    caxis(caxis_range)
    
    % colorbar set
    cb=colorbar('location','East');
    set(cb,'Position',[0.9 0.25 0.025 0.5]);
    set(cb,'YAxisLocation','right');
    set(cb,'YTick', caxis_range, 'YTickLabel', caxis_range, 'FontSize', 20, 'FontWeight', 'bold', 'FontName','Ubuntu Condensed');
    
    % title set
    til=title(figtitle,'Interpreter','none'),
    set(til,'Position',[1 80 90]);
    set(til,'FontSize',30);
    set(til,'FontWeight','Bold');
    set(til,'FontName','Ubuntu Light');
    
else
    tmp=get(gcf,'Position');
    if strcmp(BrainType,'Vert')
        bg_ords=[50, tmp(2) ,tmp(4)*1.2, tmp(4)*1.8];
        dorsal_ords=[0.05 0.24 0.8 1];
        ventral_ords=[0.05 -0.24 0.8 1];
        cb_ords=[0.9 0.25 0.025 0.5];
        title_ords=[1 -20 85];
    elseif strcmp(BrainType,'Hori')
        bg_ords=[50, tmp(2) ,tmp(4)*2.8, tmp(4)*1];
        dorsal_ords=[-0.18 0.1 0.85 0.85];
        ventral_ords=[0.28 0.1 0.85 0.85];
        cb_ords=[0.95 0.18 0.025 0.7];
        title_ords=[1 -120 -50];
    end
    
    % backgroud size
    set(gcf,'Position', bg_ords);
    
    % surface position and parameters
    % dorsal
    a(1)=axes('position',dorsal_ords);
    trisurf(surf.tri(tl,:),surf.coord(1,vl),surf.coord(2,vl),surf.coord(3,vl),...
        double(data(vl)),'EdgeColor','none');
    view(-90,0);
    daspect([1 1 1]); axis tight; camlight; axis vis3d off; zoom(0.95);
    lighting phong; material([0.5 0.5 0.2 0.2]); shading flat;
    caxis(caxis_range)
    
    % ventral
    a(2)=axes('position',ventral_ords);
    trisurf(surf.tri(tl,:),surf.coord(1,vl),surf.coord(2,vl),surf.coord(3,vl),...
        double(data(vl)),'EdgeColor','none');
    view(90,0);
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;zoom(0.95);
    lighting phong; material([0.5 0.5 0.2 0.2]); shading flat;
    caxis(caxis_range)
    
    % colorbar set
    cb=colorbar('location','East');
    set(cb,'Position',cb_ords);
    set(cb,'YAxisLocation','right');
    set(cb,'YTick', caxis_range, 'YTickLabel', caxis_range, 'FontSize', 10, 'FontWeight', 'bold', 'FontName','Ubuntu Condensed');
    
    % title set
    til=title(figtitle,'Interpreter','none'),
    set(til,'Position',title_ords);
    set(til,'FontSize',20);
    set(til,'FontWeight','Bold');
    set(til,'FontName','Ubuntu Light');
    

end
%     
id0=[0 0 cuv 0 0 cuv 0 0];
for i=1:length(a)
    set(a(i),'CLim',caxis_range);
    set(a(i),'Tag',['SurfStatView ' num2str(i) ' ' num2str(id0(i))]);
end

background='white';
whitebg(gcf,background);
set(gcf,'Color',background,'InvertHardcopy','off');

dcm_obj=datacursormode(gcf);
set(dcm_obj,'UpdateFcn',@SurfStatDataCursor,'DisplayStyle','window');

% set(gcf,'PaperPosition',[0.25 2.5 6 4.5]);
set(gcf,'PaperPositionMode','auto');

return

