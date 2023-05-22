function [ a, cb ] = plot_data( data, surf, title, background, mask );
% Slightly edited SurfStat script SurfStatViewData.m

% Added
materialProps = [0.3 0.3 0.06 5 1.0];

% Added
if nargin>4
    surftri1 = surf.tri(ismember(surf.tri(:,1), find(mask)) | ismember(surf.tri(:,2), find(mask)) | ...
        ismember(surf.tri(:,3), find(mask)), :);
    surftri2 = surf.tri(ismember(surf.tri(:,1), find(~mask)) | ismember(surf.tri(:,2), find(~mask)) | ...
        ismember(surf.tri(:,3), find(~mask)), :);
end

if nargin<3 
    title=inputname(1);
end
if nargin<4
    background='white';
end

% find cut between hemispheres, assuming they are concatenated
t=size(surf.tri,1);
v=size(surf.coord,2);
tmax=max(surf.tri,[],2);
tmin=min(surf.tri,[],2);
% to save time, check that the cut is half way
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

clim=[min(data),max(data)];
if clim(1)==clim(2)
    clim=clim(1)+[-1 0];
end

clf;
colormap(spectral(256));

h=0.39;
w=0.4;

r=max(surf.coord,[],2)-min(surf.coord,[],2);
w1=h/r(2)*r(1)*3/4;
h1=h/r(2)*r(1); % h/r(2)*r(3)


a(1)=axes('position',[0.055+0.05 0.62-0.03 h*3/4 w]);
if nargin>4
    trisurf(surftri1, surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
        double(data),'EdgeColor','none');
    hold on
    trisurf(surftri2, surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
        'FaceColor', [0.5 0.5 0.5], 'EdgeColor','none', 'FaceAlpha', 0.5);
    xlim([-80 -10])
else
    trisurf(surf.tri(tl,:),surf.coord(1,vl),surf.coord(2,vl),surf.coord(3,vl),...
        double(data(vl)),'EdgeColor','none');
end
view(-90,0); 
daspect([1 1 1]); axis tight; camlight; axis vis3d off;
lighting phong; material(materialProps);
lightangle(270,0)


a(2)=axes('position',[0.5-0.3-0.02 h1 w h]);
if nargin>4
    trisurf(surftri1, surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
        double(data),'EdgeColor','none');
    hold on
    trisurf(surftri2, surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
        'FaceColor', [0.5 0.5 0.5], 'EdgeColor','none', 'FaceAlpha', 0.5);
else
    trisurf(surf.tri,surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
        double(data),'EdgeColor','none');
end
view(0,90); 
daspect([1 1 1]); axis tight; camlight; axis vis3d off;
lighting phong; material(materialProps);
lightangle(0,90)

if cut<t
    a(3)=axes('position',[1-0.055-h*3/4-0.05 0.62-0.03 h*3/4 w]);
    if nargin>4
        trisurf(surftri1, surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
            double(data),'EdgeColor','none');
        hold on
        trisurf(surftri2, surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
            'FaceColor', [0.5 0.5 0.5], 'EdgeColor','none', 'FaceAlpha', 0.5);
        xlim([10 80])
    else
        trisurf(surf.tri(tr,:)-cuv,surf.coord(1,vr),surf.coord(2,vr),surf.coord(3,vr),...
            double(data(vr)),'EdgeColor','none');
    end
    view(90,0);
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;
    lighting phong; material(materialProps);
    lightangle(90,0)
else
    a(3)=axes('position',[1-0.055-h*3/4 0.62 h/r(2)*r(1)*3/4 w]);
    trisurf(surf.tri,surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
        double(data),'EdgeColor','none');
    view(180,0);
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;
    lighting phong; material(materialProps);
    lightangle(90,0)
end


a(4)=axes('position',[0.5-0.1+0.02 h1 w h]);
if nargin>4
    trisurf(surftri1, surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
        double(data),'EdgeColor','none');
    hold on
    trisurf(surftri2, surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
        'FaceColor', [0.5 0.5 0.5], 'EdgeColor','none', 'FaceAlpha', 0.5);
    zlim([-80 40])
else
    trisurf(surf.tri,surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
        double(data),'EdgeColor','none');
end
view(0,-90); 
daspect([1 1 1]); axis tight; camlight; axis vis3d off;
lighting phong; material(materialProps);
lightangle(0,-90)


a(5)=axes('position',[0.055+0.05 0.02+0.03 w1 h1]);
if nargin>4
    trisurf(surftri1, surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
        double(data),'EdgeColor','none');
    hold on
    trisurf(surftri2, surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
        'FaceColor', [0.5 0.5 0.5], 'EdgeColor','none', 'FaceAlpha', 0.5);
else
    trisurf(surf.tri,surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
        double(data),'EdgeColor','none');
end
view(180,0);
daspect([1 1 1]); axis tight; camlight; axis vis3d off;
lighting phong; material(materialProps);
lightangle(180,0)

a(6)=axes('position',[1-0.055-w1-0.05 0.02+0.03 w1 h1]);
if nargin>4
    trisurf(surftri1, surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
        double(data),'EdgeColor','none');
    hold on
    trisurf(surftri2, surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
        'FaceColor', [0.5 0.5 0.5], 'EdgeColor','none', 'FaceAlpha', 0.5);
else
    trisurf(surf.tri,surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
        double(data),'EdgeColor','none');
end
view(0,0);
daspect([1 1 1]); axis tight; camlight; axis vis3d off;
lighting phong; material(materialProps);
lightangle(0,0) 
    
id0=[0 0 cuv 0 0 cuv 0 0];
for i=1:length(a)
    set(a(i),'CLim',clim);
    set(a(i),'Tag',['SurfStatView ' num2str(i) ' ' num2str(id0(i))]);
end


cb=colorbar('location','South');
set(cb,'Position',[0.35 0.085 0.3 0.03]);
set(cb,'XAxisLocation','bottom');
h=get(cb,'Title');

whitebg(gcf,background);
set(gcf,'Color',background,'InvertHardcopy','off');

dcm_obj=datacursormode(gcf);
set(dcm_obj,'UpdateFcn',@SurfStatDataCursor,'DisplayStyle','window');

set(gcf,'PaperPosition',[0.25 2.5 6 4.5]);

return
end
