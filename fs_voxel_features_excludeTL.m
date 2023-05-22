function [output] = fs_voxel_features_excludeTL(subject, libdir, iso2meshdir, excludeRhTL, excludeLhTL, saveto, side)
%
% Input:
% - subject: char, path to FreeSurfer folder
% - libdir: char, path to lib-folder: https://github.com/cnnp-lab/CorticalFoldingAnalysisTools/tree/master/lib
% - iso2meshdir: char, path to ISO2MESH libaray: http://iso2mesh.sourceforge.net/
%   Licence: CC-BY
% - excludeRhTL: Whether to exclude right hemisphere temporal lobe, 0 or 1
% - excludeLhTL: Whether to exclude left hemisphere temporal lobe, 0 or 1
% - saveto: folder for output
% - side: 'l' or 'r', default is both
%
% Output: table for each hemisphere, contains values at each point on the
% pial surface
% - At (pial area) corrected using Gaussian curvature of convex hull
% - Ae (smooth pial area) corrected using Gaussian curvature of convex hull
% - At (pial area) raw
% - Ae (smooth pial area) raw
% - Average cortical thickness
% - Gaussian curvature
% - K
% - I
% - S

addpath(libdir)
addpath([libdir '/FSmatlab/'])
addpath(iso2meshdir)

if nargin < 7
    side = 'rl';
end

output = struct;

% Spherical radius around point (mm)
r = 30;

[~, SID] = fileparts(subject);

load('LUT_lobes.mat') % https://github.com/cnnp-lab/CorticalFoldingAnalysisTools/tree/master/Lobes/LUT_lobes.mat
LUT_lobes_TL = LUT_lobes(LUT_lobes(:,2)==3,1);

for hemisphere = 1:length(side)

    % Path to subject's files
    pathpre = [subject '/surf/' side(hemisphere)];
    outpath = [saveto '/' SID '/surf/' side(hemisphere)];

    [thickness, ~]  = read_curv([pathpre, 'h.thickness']);
    [pialv,pialf]   = freesurfer_read_surf([pathpre, 'h.pial']);
    [opialv,opialf] = freesurfer_read_surf([pathpre, 'h.pial-outer-smoothed']);
    [midv,midf] = freesurfer_read_surf([pathpre, 'h.mid']);
    [~,labelDK,colortable] = read_annotation([subject '/label/' side(hemisphere) 'h.aparc.annot']);
    insula = labelDK == colortable.table(string(colortable.struct_names) == "insula",5);
    
    output.([side(hemisphere) 'h']) = zeros(length(pialv),9);
    
    % Downsample pial surface and find nearest point on pial for each point on ds pial
    [midv_ds,pialf_ds] = meshresample(midv,midf,0.05); % pialf_ds = midf_ds
    pialv = single(pialv);
    [label_ds,pialv_ds,~,nearest_pialds]=matchSurfLabel_mid2pial(labelDK,pialv,midv,midv_ds);
    pialv_ds = single(pialv_ds);

    % Downsample smooth pial
    [opialv_ds, opialf_ds] = meshresample(opialv,opialf,0.1);
    clear opialv opialf

    TotalArea = NaN(length(pialv_ds),1);
    SmoothArea = NaN(length(pialv_ds),1);
    AvgThickness = NaN(length(pialv_ds),1);
    GaussCurv = NaN(length(pialv_ds),1);
%     GaussCurvPial = NaN(length(pialv_ds),1);
    
    %%
    % Only keep points on ds pial which are ot on the CC
    pialf_ds = pialf_ds((label_ds(pialf_ds(:,1)) ~= 0) & (label_ds(pialf_ds(:,2)) ~= 0) & ...
        (label_ds(pialf_ds(:,3)) ~= 0),:);

    % Only keep points on ds pial which are not on the TL
    if (side(hemisphere) == 'l' && excludeLhTL) || (side(hemisphere) == 'r' && excludeRhTL)
        pialf_ds = pialf_ds((~ismember(label_ds(pialf_ds(:,1)), LUT_lobes_TL)) & (~ismember(label_ds(pialf_ds(:,2)), LUT_lobes_TL)) & ...
            (~ismember(label_ds(pialf_ds(:,3)), LUT_lobes_TL)),:);
    end
    
    % Points which to use to get back to pial later and over which to iterate
    usepoints = unique(pialf_ds);

    % Find nearest point on ds mid for each point in mid
    % Use same relation to go back from ds pial to pial (and at the same time also sphere)
    midv = single(midv);
    nearest_sphere = match_surf(midv_ds, midv);
    
    %%
    % Label the opial_ds with DK
    % for each point on ds opial the point in pial that is
    % nearest
    nearest_opialds = match_surf(pialv, opialv_ds);
    
    %%
    % Only keep point on ds opial which are not on the CC
    label_ods = labelDK(nearest_opialds);
    opialf_ds = opialf_ds((label_ods(opialf_ds(:,1)) ~= 0) & (label_ods(opialf_ds(:,2)) ~= 0) & ...
        (label_ods(opialf_ds(:,3)) ~= 0),:);
    
    % Only keep points on ds opial which are not on the TL
    if (side(hemisphere) == 'l' && excludeLhTL) || (side(hemisphere) == 'r' && excludeRhTL)
        opialf_ds = opialf_ds((~ismember(label_ods(opialf_ds(:,1)), LUT_lobes_TL)) & (~ismember(label_ods(opialf_ds(:,2)), LUT_lobes_TL)) & ...
            (~ismember(label_ods(opialf_ds(:,3)), LUT_lobes_TL)),:);
    end
    
    % Nearest point on ds pial for points in ds smooth pial
    nearest_opialds = match_surf(pialv_ds, opialv_ds);

    %%
    % Find thickness of vertices of pial closest to each ds vertex
    thickness_ds = thickness(nearest_pialds);
    thickness_ds_smoothed = NaN(length(pialv_ds),1);

    for index = 1:length(usepoints)
        point = usepoints(index);
        neighbourIDs = point;
        vertex = pialv_ds(neighbourIDs,:);
        rmin = 0;
        while rmin < r
            % Find neighbours of neighbours
            fid = ismember(pialf_ds(:,1),neighbourIDs) | ...
                ismember(pialf_ds(:,2),neighbourIDs) | ...
                ismember(pialf_ds(:,3),neighbourIDs);

            new = setdiff(unique(pialf_ds(fid,:)), neighbourIDs);
            d = sqrt(sum((pialv_ds(new,:)-vertex).^2,2));
            rmin = min(d);
            new = new(:);
            neighbourIDs = [neighbourIDs; new];
        end
        thickness_ds_smoothed(point) = mean(thickness_ds(neighbourIDs), 'omitnan');
    end

    % Thickness of faces in ds pial
    ThicknessFB = makeFacebased(thickness_ds, pialf_ds); % Don't use averaged version since the faces are being averaged later
    
    GruPial = getGaussianCurvPart(pialf_ds,pialv_ds,1:length(pialv_ds));

    %%
    for index = 1:length(usepoints)

        point = usepoints(index);

        % Check if it is on the midline or if the thickness is too small
        % (for healthy adults)
        if thickness_ds_smoothed(point) < 1.8 || label_ds(point) == 0
            continue;
        else

            % Index of point being looked at
            neighbourIDs = point;

            % Coordinates of that point
            vertex = pialv_ds(neighbourIDs,:);

            % Repeat until new neighbours are all > r mm from the point
            rmin = 0;
            while rmin < r

                % Find neighbours of neighbours
                fid = ismember(pialf_ds(:,1),neighbourIDs) | ...
                    ismember(pialf_ds(:,2),neighbourIDs) | ...
                    ismember(pialf_ds(:,3),neighbourIDs);

                % Check how far all new ones are from vertex, add the ones that are
                % close enough
                new = setdiff(unique(pialf_ds(fid,:)), neighbourIDs);
                d = sqrt(sum((pialv_ds(new,:)-vertex).^2,2));
                rmin = min(d);
                new = new(:);

                neighbourIDs = [neighbourIDs; new(d<r)];
            end
            
            % Label ds pial
            label = zeros(length(pialv_ds),1);
            label(neighbourIDs) = 1;
            
            % Find edges of patch
            isBoundary = zeros(length(neighbourIDs),1);
            for kl = 1:length(neighbourIDs)
                fid = pialf_ds(:,1)==neighbourIDs(kl) | ...
                      pialf_ds(:,2)==neighbourIDs(kl) | ...
                      pialf_ds(:,3)==neighbourIDs(kl);

                v_neigh = unique(pialf_ds(fid,:));
                ncolors = numel(unique(label(v_neigh)));
                if ncolors > 1
                    isBoundary(kl) = 1;
                end
            end
            edge = neighbourIDs(isBoundary==1);
            
            %%
            % Fill potential holes
            if ~isempty(edge)
                newPoints = edge(1);
                points_left = 1;
                while points_left
                    % Check we don't start with a vertex that is between a hole and the
                    % edge
                    on_edge = ismember(newPoints(1),edge);
                    while on_edge
                        bndNeigh =  ismember(pialf_ds(:,1),newPoints(1)) | ...
                            ismember(pialf_ds(:,2),newPoints(1)) | ...
                            ismember(pialf_ds(:,3),newPoints(1));

                        bndNeigh = unique(pialf_ds(bndNeigh,:));
                        bndNeigh = bndNeigh(ismember(bndNeigh, neighbourIDs(isBoundary == 1)));
                        if length(bndNeigh) > 3
                            edge(1) = [];
                            if isempty(edge)
                               points_left = 0;
                               on_edge = 0;
                            else
                               newPoints = edge(1);
                            end
                        else
                            on_edge = 0;
                        end
                    end

                    % Find neighbours of neighbours
                    fid = ismember(pialf_ds(:,1),newPoints) | ...
                        ismember(pialf_ds(:,2),newPoints) | ...
                        ismember(pialf_ds(:,3),newPoints);

                    new = unique(pialf_ds(fid,:));
                    edge(ismember(edge,new)) = [];
                    new = setdiff(new, [neighbourIDs; newPoints]);
                    new = new(:);
                    newPoints = [newPoints; new];
                    newPoints = setdiff(newPoints, neighbourIDs);
                    d = sqrt(sum((pialv_ds(newPoints,:)-vertex).^2,2));

                    if isempty(new)
                        neighbourIDs = [neighbourIDs; newPoints];
                        isBoundary = [isBoundary; zeros(length(newPoints),1)];
                        if isempty(edge)
                           points_left = 0;
                        else
                           newPoints = edge(1);
                        end
                    elseif max(d) > 50
                        % Run only on edge until no new points
                        newEdge = 1;
                        while ~isempty(newEdge)
                            fid = ismember(pialf_ds(:,1),newPoints) | ...
                                ismember(pialf_ds(:,2),newPoints) | ...
                                ismember(pialf_ds(:,3),newPoints);

                            new = unique(pialf_ds(fid,:));

                            newEdge = new(ismember(new, edge));
                            edge(ismember(edge,new)) = [];

                            new = setdiff(new, [neighbourIDs; newPoints]);
                            new = new(:);
                            newPoints = [newPoints; new];
                        end
                        if isempty(edge)
                           points_left = 0;
                        else
                           newPoints = edge(1);
                        end
                    end
                end
            end

            %%
            % Relabel ds pial
            label = zeros(length(pialv_ds),1);
            label(neighbourIDs) = 1;

            % Total area
            % Only use area of faces that are completely within circle
            [liaa] = ismember(pialf_ds, neighbourIDs);
            aid = (sum(liaa,2) == 3);

            TotalAreai = zeros(length(aid), 1);
            TotalAreai(aid) = calcTriangleArea(pialf_ds(aid,:), pialv_ds);
            TotalArea(point) = sum(TotalAreai);

            % Average thickness in circle
            AvgThickness(point) = sum(ThicknessFB(aid>0).*TotalAreai(aid>0))/TotalArea(point);
            
            % Label the ds smooth pial
            label_smooth_ds = label(nearest_opialds);
            
            % Find the smooth area
            sids = find(label_smooth_ds == 1);
            [liaa] = ismember(opialf_ds,sids);
            aid = (sum(liaa,2) == 3);

            SmoothArea(point) = sum(calcTriangleArea(opialf_ds(aid,:),opialv_ds));
 
            % Gauss. curv. of patch
            ov_ids = find(label_smooth_ds == 1);
            Lobepoints = opialv_ds(ov_ids,:);

            if length(Lobepoints) > 4

                CHSf_ds = convhull(Lobepoints);

                Gru = getGaussianCurvPart(CHSf_ds,Lobepoints,1:length(ov_ids));

                isBoundary = zeros(length(ov_ids),1);

                for kl = 1:length(ov_ids)
                    fid = opialf_ds(:,1)==ov_ids(kl) | ...
                          opialf_ds(:,2)==ov_ids(kl) | ...
                          opialf_ds(:,3)==ov_ids(kl);

                    v_neigh = unique(opialf_ds(fid,:));
                    ncolors = numel(unique(label_smooth_ds(v_neigh)));
                    if ncolors > 1
                        isBoundary(kl) = 1;
                    end
                end

                GaussCurv(point)  = sum(Gru(isBoundary == 0));
            end


%             % Gauss. curv. of patch from pial
%             Lobepoints = pialv_ds(neighbourIDs,:);
% 
%             if length(Lobepoints) > 4
% 
%                 GruNeigh = GruPial(neighbourIDs);
% 
%                 isBoundary = zeros(length(neighbourIDs),1);
% 
%                 for kl = 1:length(neighbourIDs)
%                     fid = pialf_ds(:,1)==neighbourIDs(kl) | ...
%                           pialf_ds(:,2)==neighbourIDs(kl) | ...
%                           pialf_ds(:,3)==neighbourIDs(kl);
% 
%                     v_neigh = unique(pialf_ds(fid,:));
%                     ncolors = numel(unique(label(v_neigh)));
%                     if ncolors > 1
%                         isBoundary(kl) = 1;
%                     end
%                 end
% 
%                 GaussCurvPial(point)  = sum(GruNeigh(isBoundary == 0));
%             end
        end
    end

    %%
    % Correct surface areas
    At_dash = TotalArea.*4*pi./GaussCurv;
    Ae_dash = SmoothArea.*4*pi./GaussCurv;
    
    % Set values where CH GC is below 0.16 to NaN
    At_dash(GaussCurv < 0.16) = NaN;
    Ae_dash(GaussCurv < 0.16) = NaN;
    
    % Get volume as T*At (uncorrected)
    volume = TotalArea.*AvgThickness;


    %% Data imputation
    % Find points that are NaN and not excluded intentionally
    missingValues = find(isnan(At_dash) | isnan(Ae_dash) | At_dash == Inf | Ae_dash == Inf | At_dash == -Inf | Ae_dash == -Inf);
    missingValues = missingValues(ismember(missingValues, usepoints));
    
    At_dashTMP = NaN(size(At_dash));
    Ae_dashTMP = NaN(size(Ae_dash));

    % Loop over missing vertices and take avg of surrounding values (but
    % put into tmp so I don't reuse the same values again)
    for kl = 1:length(missingValues)
        fid = pialf_ds(:,1)==missingValues(kl) | ...
              pialf_ds(:,2)==missingValues(kl) | ...
              pialf_ds(:,3)==missingValues(kl);

        v_neigh = unique(pialf_ds(fid,:));
        
        % Make sure I don't overwrite anything if not all measures are
        % missing
        if isnan(At_dash(missingValues(kl))) || At_dash(missingValues(kl)) == Inf || At_dash(missingValues(kl)) == -Inf
            At_dashTMP(missingValues(kl)) = mean(At_dash(v_neigh), 'omitnan');
        end
        if isnan(Ae_dash(missingValues(kl))) || Ae_dash(missingValues(kl)) == Inf || Ae_dash(missingValues(kl)) == -Inf
            Ae_dashTMP(missingValues(kl)) = mean(Ae_dash(v_neigh), 'omitnan');
        end
    end

    At_dash(missingValues) = At_dashTMP(missingValues);
    Ae_dash(missingValues) = Ae_dashTMP(missingValues);

    K = log10(At_dash) - 5/4*log10(Ae_dash) + ...
        1/2*log10(AvgThickness);
    I = log10(At_dash) + log10(Ae_dash) + ...
        2*log10(AvgThickness);
    S = 3/2*log10(At_dash) + 3/4*log10(Ae_dash) - ...
        2*9/4*log10(AvgThickness);

    % Go back to pial surface
    volume = volume(nearest_sphere);

    output.([side(hemisphere) 'h'])(:,1) = At_dash(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,2) = Ae_dash(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,3) = TotalArea(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,4) = SmoothArea(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,5) = AvgThickness(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,6) = GaussCurv(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,7) = K(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,8) = I(nearest_sphere);
    output.([side(hemisphere) 'h'])(:,9) = S(nearest_sphere);
    
    % Remove Inf & insula
    output.([side(hemisphere) 'h'])(output.([side(hemisphere) 'h']) == Inf | output.([side(hemisphere) 'h']) == -Inf) = NaN;
    output.([side(hemisphere) 'h'])(insula,:) = NaN;

    fnum = length(pialf);
    
    if (side(hemisphere) == 'l' && excludeLhTL) || (side(hemisphere) == 'r' && excludeRhTL)
        write_curv([outpath 'h.PialAreaNoTL'], output.([side(hemisphere) 'h'])(:,1), fnum);
        write_curv([outpath 'h.SmoothPialAreaNoTL'], output.([side(hemisphere) 'h'])(:,2), fnum);
        write_curv([outpath 'h.AvgCortThicknessNoTL'], output.([side(hemisphere) 'h'])(:,5), fnum);
        write_curv([outpath 'h.VolumeNoTL'], volume, fnum);
    else
        write_curv([outpath 'h.PialArea'], output.([side(hemisphere) 'h'])(:,1), fnum);
        write_curv([outpath 'h.SmoothPialArea'], output.([side(hemisphere) 'h'])(:,2), fnum);
        write_curv([outpath 'h.AvgCortThickness'], output.([side(hemisphere) 'h'])(:,5), fnum);
        write_curv([outpath 'h.Volume'], volume, fnum);
    end
   
    output.([side(hemisphere) 'h']) = array2table(output.([side(hemisphere) 'h']), 'VariableNames',{'PialArea','SmoothPialArea','PialAreaRaw','SmoothPialAreaRaw', ...
        'AvgCortThickness','GaussCurv','K','I','S'});
end

% Write to file
if excludeLhTL && excludeRhTL
    save([saveto '/' SID '/surf/voxel_features_excludedBothTL.mat'], 'output')
elseif excludeLhTL
    save([saveto '/' SID '/surf/voxel_features_excludedLhTL.mat'], 'output')
elseif excludeRhTL
    save([saveto '/' SID '/surf/voxel_features_excludedRhTL.mat'], 'output')
else
    save([saveto '/' SID '/surf/voxel_features.mat'], 'output')
end
end
