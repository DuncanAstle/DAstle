% GHSOM_DATALABES  Assign the labels of the dataitems to their best matching unit.
%
%  ghMap = ghsom_datalabels(ghMap, sData, [mode])
%
%  ghMap = ghsom_train(sData);
%  ghMap = ghsom_datalabes(ghMap, sData);
%  ghVisu = ghsom_visualize_grid(ghMap,'layer');
%  ghMap = ghsom_visualize_labels(ghMap, ghVisu);
%
%  mode (string) can be omitted when set to 'vote' only the most frequent label
%  is used, see som_autolabel for details.
%
% See also GHSOM_VISUALIZE_LABELS, GHSOM_VISUALIZE_GRID

% Copyright (c) 2001 by Elias Pampalk

% Version 1.0 Elias Pampalk 17072001

% Functions used of the SOM Toolbox:
%  SOM_LABEL
%  SOM_AUTOLABEL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ghMap = ghsom_datalabels(ghMap, sData, mode)

% label the map
if nargin==2,
    for i=1:length(ghMap.parent),
        ghMap.sMap{i}=som_label(ghMap.sMap{i}, 'clear', 'all');
        ghMap.sMap{i}=som_autolabel(ghMap.sMap{i}, ghsom_idx2data(sData,ghMap.dataitems{i}));
    end
else
    for i=1:length(ghMap.parent),
        ghMap.sMap{i}=som_label(ghMap.sMap{i}, 'clear', 'all');
        ghMap.sMap{i}=som_autolabel(ghMap.sMap{i}, ghsom_idx2data(sData,ghMap.dataitems{i}),mode);
    end
end