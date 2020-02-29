% GHSOM_READ_CODS  Read GHSOM from SOM_PAK format.
%
%  ghMap = ghsom_read_cods(directory)
%  ghMap = ghsom_map_data(ghMap, sData)
%
% Note: SOM_READ_COD is used.
%
% See also: GHSOM_WRITE_CODS, GHSOM_MAP_DATA

% Copyright (c) 2001 by Elias Pampalk

% Version 1.0 Elias Pampalk 17072001

function ghMap = ghsom_read_cods(directory)

if ~strcmp(directory(end),'/'),
    directory = [directory, '/'];
end

d=dir(directory);

ghMap.expand_units=[]; % dont bother to recalculate these
ghMap.expand_units_qe=[];
ghMap.dataitems=[];

% read files...
for i=1:length(d),
    if ~d(i).isdir,
        if strcmp(d(i).name(end-3:end),'.cod'), % only take files with extension '.cod'
            [str, R] = strtok(d(i).name(1:end-4),'_');
            midx = str2num(str);
            [str, R] = strtok(R,'_'); 
            ghMap.layer(midx) = str2num(str);
            [str, R] = strtok(R,'_');
            ghMap.parent(midx) = str2num(str);
            [str, R] = strtok(R,'_');
            ghMap.parent_unit(midx) = str2num(str);
            ghMap.sMap{midx}=som_read_cod([directory,d(i).name]);
        end
    end
end
