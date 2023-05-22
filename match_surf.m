function [midx] = match_surf(labelled_surf, unlabelled_surf)

steps=1000;

np=size(unlabelled_surf,1);
midx=NaN(1,size(unlabelled_surf,1));

for i=1:steps:np
    if i+steps-1>np
        ei=np;
    else
        ei=i+steps-1;
    end
    
    cp=unlabelled_surf(i:ei,:);
    
    d1=labelled_surf(:,1)-cp(:,1)';
    d2=labelled_surf(:,2)-cp(:,2)';
    d3=labelled_surf(:,3)-cp(:,3)';
    
    d=sqrt(d1.^2+d2.^2+d3.^2);
    
    [~,midx(i:ei)]=min(d);
end
end