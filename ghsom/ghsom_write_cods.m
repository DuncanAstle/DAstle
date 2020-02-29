% GHSOM_WRITE_CODS  Store GHSOM in SOM_PAK format.
%
%  ghsom_write_cods(ghMap, directory)
%
% Note: SOM_WRITE_COD is used.
%
% The files created are of the SOM_PAK format and contain the
% codebooks for one sMap.
% The file names are defined as:
% directory + map idx + layer + parentmap + unit off the parent map + '.cod'
%
% e.g.: 
% './1_1_0_1.cod' first layer
% './2_2_1_2.cod' 2nd layer, expanded from 2nd unit on 1st map
% './3_2_1_3.cod' 2nd layer, expanded from 3rd unit on 1st map
% './4_3_3_4.cod' 3rd layer, expanded from 4th unit on 3rd map

% Copyright (c) 2001 Elias Pampalk

% Version 1.0 Elias Pampalk 17072001

function ghsom_write_cods(ghMap,directory)

if ~strcmp(directory(end),'/'),
    directory = [directory, '/'];
end

for i=1:length(ghMap.parent),
    filename=[directory,...
            num2str(i),'_',...
            num2str(ghMap.layer(i)),'_',...
            num2str(ghMap.parent(i)),'_',...
            num2str(ghMap.parent_unit(i)),...
            '.xls'];
    % because of som toolbox bug remove labels from map!
    ghMap.sMap{i}=som_label(ghMap.sMap{i}, 'clear', 'all');
    
    som_write_cod(ghMap.sMap{i},filename);
end         
