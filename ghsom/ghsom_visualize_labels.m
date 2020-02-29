% GHSOM_VISUALIZE_LABELS  Plot labels onto grid.
%
%  ghsom_visualize_labels(ghMap, ghVisu)
%
% Example using the datalabls:
%  ghMap = ghsom_train(Data);
%  ghMap = ghsom_datalabels(ghMap);
%  ghVisu = ghsom_visualize_grid(ghMap,'layer');
%  ghMap = ghsom_visualize_labels(ghMap, ghVisu);
%
% Also for use with WebSOM and LabelSOM labelling.
%
% See also GHSOM_DATALABELS, GHSOM_LABELSOM, GHSOM_WEBSOM_LABEL

% Copyright (c) 2001 by Elias Pampalk

% Version 1.0 Elias Pampalk 17072001

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ghsom_visualize_labels(ghMap, ghVisu)

for map_idx=1:length(ghMap.sMap),
    for unit_idx=1:prod(ghMap.sMap{map_idx}.topol.msize),
        if no_sub_unit(map_idx,unit_idx,ghMap.parent,ghMap.parent_unit), 
            text(...
                ghVisu.coordinates_xy{map_idx}(unit_idx,1),...
                ghVisu.coordinates_xy{map_idx}(unit_idx,2),...
                ghMap.sMap{map_idx}.labels(unit_idx,...
                find(strcmp(ghMap.sMap{map_idx}.labels(unit_idx,:),'')==0)),...
                'HorizontalAlignment','center','fontsize',6);
        end
    end    
end

function is_true=no_sub_unit(map_idx, unit_idx, parentmap, parent_unit)
is_true=1;
children=find(parentmap==map_idx);
if ~isempty(children),
    if ~isempty(find(parent_unit(children)==unit_idx)),
        is_true=0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%