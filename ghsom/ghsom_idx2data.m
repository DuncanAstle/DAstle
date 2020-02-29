% GHSOM_IDX2DATA  Having the complete sData set, get the sub set which belongs to a submap. 
%
%  subsData = ghsom_idx2data(sData, dataitems)
%
%  'dataitems' (vector) which dataitems belong to the sub sData set. (ghMap.dataitems)
%  'sData'     (struct) the complete data set. (SOM Toolbox format)
%  'subsData'  (struct) the sub set. (SOM Toolbox format)
%
% See also GHSOM_TRAIN, GHSOM_VISUALIZE_DATALABELS

% Copyright (c) 2001 by Elias Pampalk

% Version 1.0 Elias Pampalk 17072001


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sSubData = ghsom_idx2data(sData, dataitems)
% get an sData structure knowing which subset (dataitems) from the
% complete sData are wanted.

sSubData.type        = 'som_data';
sSubData.data        = sData.data(dataitems,:);
sSubData.labels      = sData.labels(dataitems,:);
sSubData.name        = [sData.name,' (subset)'];
sSubData.comp_names  = sData.comp_names;
sSubData.comp_norm   = sData.comp_norm;
sSubData.label_names = []; % ?
%% END Function: find_expand_units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
