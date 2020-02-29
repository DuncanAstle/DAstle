% GHSOM_LABELSOM  Adds labels to the GHSOM based on the LabelSOM method.
%
%  ghMap = ghsom_labelsom(ghMap, sData, threshold, mlabels)
%
% Input and Output arguments: 
%   ghMap      (struct) ghsom Toolbox structure for a GSHOM map.
%   sData      (struct) SOM Toolbox's structure for a SOM formatted Data.
%   threshold  (scaler) minimum value in the codebook vector 
%                       element, (to remove elements representing 0 or null 
%                       features).
%   mlabels (scaler) maximum number of labels to be placed on a unit.
%   
% ghMap = ghsom_labelsom(ghMap,sData,0.2,5)
%
% FUNCTIONS USED:
%              som_bmus
%
% NOTE:        Typical value for threshold is 0.1 or 0.2 (For removing 
%              features that have a zero value which consitutes to null 
%              feature in the vector [Only applicable to text info]). 
%              This is typical for normalised codebook vectors between
%              0 to 1. 

% Copyright (c) 2001 by Elias Pampalk, Alvin Chan

% Version 1.0 Alvin Chan 23102000 30042001
% Version 2.0 Elias Pampalk 17072001

function ghMap = ghsom_labelsom(ghMap, sData, threshold, mlabels)

for m=1:length(ghMap.parent),
    
    sMap = ghMap.sMap{m};
    subData = ghsom_idx2data(sData,ghMap.dataitems{m});
    
    [munits dims] = size(sMap.codebook);
    
    % 0. remove old labels
    sMap.labels=repmat({''},munits,mlabels);
    
    % 1. calculate the mean quantisation error for each map unit
    bmus = som_bmus(sMap, subData);
    unit_std = zeros(munits,dims);
    unit_mapped = zeros(munits,1);
    
    for i = 1:munits,
        if any(bmus==i),
            unit_std(i,:) = std(subData.data(bmus==i,:));
        end
        unit_mapped(i) = sum(bmus==i);
    end
    
    % 2. add labels
    for i = 1:munits,
        if unit_mapped(i),
            idx1 = find(sMap.codebook(i,:) >= threshold);
            [dummy, idx2] = sort(unit_std(i,idx1,:)); % want smalles std's
            
            idx12=idx1(idx2);
            if ~isempty(idx12),            
                sMap.labels(i,1:min(mlabels,length(idx12))) = subData.comp_names(idx12(1:min(mlabels,length(idx12))));
            end
        end
    end
    ghMap.sMap{m}=sMap;
end
