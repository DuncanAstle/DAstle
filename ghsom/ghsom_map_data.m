% GHSOM_MAP_DATA  Map sData onto a ghMap (set ghMap.dataitems).
%
%  ghMap = ghsom_map_data(ghMap, sData)
%
% See also: SOM_DIVIDE

% Copyright (c) 2001 by Elias Pampalk

% Version 1.0 Elias Pampalk 17072001

function ghMap = ghsom_map_data(ghMap, sData)

ghMap.dataitems{1}=1:(size(sData.data,1));

for i=2:length(ghMap.parent),
    p = ghMap.parent(i);
    [dummy ghMap.dataitems{i}] = som_divide( ...
        ghMap.sMap{p}, ...
        ghsom_idx2data(sData,ghMap.dataitems{p}), ... 
        ghMap.expand_units{p}(i), ...
        'index');
    ghMap.dataitems{i} = ghMap.dataitems{p}(ghMap.dataitems{i});
end