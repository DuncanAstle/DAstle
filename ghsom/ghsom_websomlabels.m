% GHSOM_WEBSOMLABELS  Adds labels to a map based on Lagus and Kaskis method.
%
%  ghMap = ghsom_websomlabels(ghMap, sData, max_labels)
%
% Input and Output arguments: 
%   ghMap      (struct) ax_ghsom Toolbox structure for a GSHOM map.
%   sData      (struct) SOM Toolbox's structure for a SOM formatted Data.
%   max_labels (scaler) maximum number of labels to be placed on a unit.
%   
% ghMap = ghsom_websomlabels(ghMap,sData,5)

% Copyright (c) 2001 Elias Pampalk, Alvin Chan

% Version 1.0 Alvin Chan 23102000 30042001
% Version 2.0 Elias Pampalk 17072001

function ghMap = ghsom_websomlabels(ghMap, sData, max_labels)
for m=1:length(ghMap.parent),
    
    sMap = ghMap.sMap{m};
    subData = ghsom_idx2data(sData,ghMap.dataitems{m});
    
    [munits dims] = size(sMap.codebook);
    
    % 0. remove old labels
    sMap.labels=repmat({''},munits,max_labels);
    
    % 1. normalize codebook values by words in codebook (F_cluster)
    cluster=repmat(sum(sMap.codebook,2),1,dims);
    cluster(cluster==0)=eps;
    cb = sMap.codebook./cluster;
    
    % 2. calculate G_0 (normalize by frequency of word in collection)
    collection=repmat(sum(cb,1),munits,1);
    collection(collection==0)=eps;
    G_0 = cb.^2./collection;
    
    % 2. add labels
    for i = 1:munits,
        [dummy, idx] = sort(G_0(i,:)); % want biggest values
        sMap.labels(i,:) = subData.comp_names(idx(end-max_labels+1:end));
    end
    ghMap.sMap{m}=sMap;
end