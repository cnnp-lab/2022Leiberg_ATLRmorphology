
%%% Requires surfstat library https://www.math.mcgill.ca/keith/surfstat/

addpath('surfstat')

clear all

% folder with patient data that has been age and sex corrected (using gam),
% hemisphere of right onset patients have been flipped, 
% left and right hemispheres of all patients have been collated
inputFolder = 'Data/age_sex_corrected/';

% save outputs and figs
outputFolder = 'figures/';

% fs avg surf to plot
surf = SurfStatReadSurf({'fsaverage/surf/lh.pial', 'fsaverage/surf/rh.pial'});

measures = ["T" "At" "Ae" "K" "I" "S"];
KIS = importdata('KIS_loadings_AtAeTxKIS.csv');

pThresh = 0.05;

%%
for vari = 1:6

    outputs = struct;

    if vari < 4
        load([inputFolder char(measures(vari)) '_allPatients.mat'])
    
        dataPre = allPatients.dataPre;
        dataPost = allPatients.dataPost;
    else
        load([inputFolder char(measures(1)) '_allPatients.mat'])
        TdataPre = allPatients.dataPre;
        TdataPost = allPatients.dataPost;
        
        load([inputFolder char(measures(2)) '_allPatients.mat'])
        AtdataPre = allPatients.dataPre;
        AtdataPost = allPatients.dataPost;
        
        load([inputFolder char(measures(3)) '_allPatients.mat'])
        AedataPre = allPatients.dataPre;
        AedataPost = allPatients.dataPost;

        % Compute KIS
        dataPre = AtdataPre*KIS(1,vari-3) + AedataPre*KIS(2,vari-3) + TdataPre*KIS(3,vari-3);
        dataPost = AtdataPost*KIS(1,vari-3) + AedataPost*KIS(2,vari-3) + TdataPost*KIS(3,vari-3);
    end

    % Centre data subject-wise for paired analysis
    [dataPre, dataPost] = centre_data(dataPre, dataPost);

    subN = size(dataPre,1);

    % Drop vertices where there is missing data
    dataPre(:,(sum(isnan(dataPre))>0) | (sum(isnan(dataPost))>0)) = NaN;
    dataPost(:,(sum(isnan(dataPre))>0) | (sum(isnan(dataPost))>0)) = NaN;

    data = [dataPre; dataPost];

%%
    pvalsPos = NaN(1, size(surf.coord,2));
    pvalsNeg = NaN(1, size(surf.coord,2));
    
    % Set up contrast and design matrix
    grouping = [repmat({'pre'}, subN,1); repmat({'post'}, subN,1)];
    grouping = term(grouping);
    contrast = grouping.post - grouping.pre;
    M = grouping;

    % Run surfstat
    slm = SurfStatLinMod( data, M, surf );
    slm = SurfStatT( slm, contrast );
    pval = SurfStatP( slm );

    % F-squared
    Rab = 1 - slm.SSE./sum((data - mean(data)).^2);
    effects = Rab(1,:)./(1-Rab(1,:));
    
    % Cluster p-values
    if isfield(pval, 'C') % It only computes cluster p-values if there are t-values above treshold
        pvalsPos = pval.C;
    end

    % Contrast other way around to test for negative effects
    contrast = -contrast;
    
    slm = SurfStatT( slm, contrast );        
    pval = SurfStatP( slm );
    
    % Cluster p-values
    if isfield(pval, 'C') % It only computes cluster p-values if there are t-values above treshold
        pvalsNeg = pval.C;
    end
 
%% Plots
    if sum(isnan(pvalsPos),'all') + sum(isnan(pvalsNeg),'all') ~= 2*size(surf.coord,2)

        % Threshold p-values
        signifclus=pvalsPos<pThresh;        
        signifclus2=pvalsNeg<pThresh;
        signifclusBoth = logical(signifclus+signifclus2);

        mask = ~isnan(data(1,:)); % Regions that were excluded s.a. TL

        outputs.significlus = signifclusBoth;
        outputs.effects = effects;
        outputs.pvalsPos = pvalsPos;
        outputs.pvalsNeg = pvalsNeg;
        outputs.mask = mask;
        outputs.dataPre = dataPre;
        outputs.dataPost = dataPost;

        effects(~signifclusBoth) = NaN;
        effects(signifclus2) = -effects(signifclus2);

        save([outputFolder '/outputs_' char(measures(vari)) '_LTLEandRTLE.mat'], 'outputs');

        if sum(signifclusBoth) ~= 0

            % Cap at 0.5 for plots
            effects(effects > 0.5) = 0.5;
            effects(effects < -0.5) = -0.5;

            plot_effect(effects, signifclusBoth, mask, char(measures(vari)), surf, outputFolder);
    
            close all
        end
    end
end