% GHSOM_VISUALIZE_GRID  Plot a 2D grid representing the structure of the GHSOM. 
%
%  ghVisu = ghsom_visualize_grid(ghMap, [color, [component_planes or #layers]])
%
%  ghMap = ghsom_train(Data);
%  ghVisu = ghsom_visualize_grid(ghMap);
%
%  color (string)
%                   'none'      (default), no coloring is done.
%                   'layer'     the deepest layer is white, the top layers are gray
%                               the idea is that clusters are areas with high data
%                               densities on the map. High density areas are hierarchically
%                               expanded by the SOM. (Original idea by Alvin Chan)
%                               Use layer to determine how many layers of the map shall be
%                               visualized (default = 0 all layers);
%                   'component' visualize component planes. 
%  component_planes (vector)  which component planes to visualize. (subplots are used)
%
%  Returns a structure containing information on the coordinates ('coordinates_xy')
%  and the scaleing factors ('unit_size_xy'). Both are needed to add for example
%  labels to the grid.
%
% See also GHSOM_TRAIN

% Copyright (c) 2001-2002 by Elias Pampalk, Alvin Chan

% Version 1.0 Elias Pampalk 17072001
% Version 1.1 Elias Pampalk 20052002
% Version 1.2 Elias Pampalk 31052002

function v = ghsom_visualize_grid(ghMap, color, component_planes)
% Note: msize defined as [rows columns] X = columne, Y = row

max_layer = ghMap.layer(end); % should be max 

if nargin==1,
    color='none';
end

if nargin==2 & strcmp(color,'component'),
    error('Error (ghsom_visualize): Which component(s)?');
end

num_layers = 0;
plots=1;
if strcmp(color,'component'),
    plots=length(component_planes);
elseif nargin==3,
    num_layers = component_planes;
end

if plots<=3,
    subplots_cols=plots; 
    subplots_rows=1;
else
    subplots_rows=ceil(sqrt(plots));
    subplots_cols=subplots_rows;
end

switch color
case 'component',
    layer_typ=0;
case 'layer',
    layer_typ=1;    
otherwise
    error('(ghsom_visualize_grid) Unknown color type.');
end

for p=1:plots,
    s=subplot(subplots_rows,subplots_cols,p); axis off; hold on;
    colors=NaN;
    if ~layer_typ
        colormap gray;
        h=title(ghMap.sMap{1}.comp_names{component_planes(p)});
        set(h,'verticalalignment','middle')
        colors=ghMap.sMap{1}.codebook(:,component_planes(p))';
    else
        mycolors=(repmat(linspace(0.7,1,64)',1,3));
        colors=mycolors(1,:); % first layer
    end
    
    % first layer (special case because no parent unit exists)
    [v.coordinates_xy{1}, v.unit_size_xy{1}] = add_map(colors, ghMap.sMap{1}.topol.msize, [], [], layer_typ);
    
    % all other layers
    for map_idx=2:length(ghMap.sMap),
        if num_layers & ghMap.layer(map_idx)>num_layers,
            break;
        end
        
        if ~layer_typ
            colors=ghMap.sMap{map_idx}.codebook(:,component_planes(p))';
        else
            colors=mycolors(round(ghMap.layer(map_idx)/max_layer*64),:);
        end
        
        [v.coordinates_xy{map_idx}, v.unit_size_xy{map_idx}] = add_map( ...
            colors, ...
            ghMap.sMap{map_idx}.topol.msize, ...                                        
            v.unit_size_xy{ghMap.parent(map_idx)}, ...                               
            v.coordinates_xy{ghMap.parent(map_idx)}(ghMap.parent_unit(map_idx),:), layer_typ);  
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coordinates_xy, unit_size_xy] = ...
    add_map(colors, msize_yx, size_parent_unit_xy, parent_coordinates_xy,layer_typ)
% WHERE:
%   msize_xy                -> map size [x, y] of current map to draw (rectangular shape)
%   size_parent_unit_xy     -> used as scaling factors [x, y]
%   parent_coordinates_xy   -> center of the parent unit [x, y] coordinates

msize_xy([2 1]) = msize_yx;

% make submaps a little smaller to recognize which maps are submaps from which maps
if layer_typ,
    SHRINKER = .9; 
else
    SHRINKER = 1;
end

% simplify usage for first layer map
if isempty(size_parent_unit_xy), 
    size_parent_unit_xy = [1,1]; 
    parent_coordinates_xy = [0,0]; 
end

% basic rectangle for patch: size 1 by 1, centered at 0,0
rectangle = [-0.5000, -0.5000; -0.5000, 0.5000; 0.5000, 0.5000; 0.5000, -0.5000]; 

% width and height of map units depend on how many units need to fit in the space suplied by parent unit
unit_size_xy = size_parent_unit_xy./msize_xy*SHRINKER;

% define mapping of the units (1..max) into coordinate space (x,y)
%   unit 1 is (1,max(y)) upper left 
%   unit 2 is (1,max(y)-1) upper left, down one
%   unit max is (max(x),1) lower right
coordinates_xy = [...
        reshape(repmat(1:msize_xy(1),msize_xy(2),1),1,prod(msize_xy)); ...
        repmat(msize_xy(2):-1:1,1,msize_xy(1))]';

% translate (move) into center (0, 0)
coordinates_xy(:,1) = coordinates_xy(:,1) - mean(coordinates_xy(:,1));
coordinates_xy(:,2) = coordinates_xy(:,2) - mean(coordinates_xy(:,2));

% scale to right size and translate (move) to right position
coordinates_xy(:,1) = coordinates_xy(:,1)*unit_size_xy(1) + parent_coordinates_xy(1);
coordinates_xy(:,2) = coordinates_xy(:,2)*unit_size_xy(2) + parent_coordinates_xy(2);

% create rectangle with appropriate scaling, centered at the coordinates calculated
rectangles_x = ...
    repmat(rectangle(:,1),1,prod(msize_xy)) * unit_size_xy(1) + ...
    repmat(coordinates_xy(:,1)',4,1);
rectangles_y = ...
    repmat(rectangle(:,2),1,prod(msize_xy)) * unit_size_xy(2) + ...
    repmat(coordinates_xy(:,2)',4,1);

patch(rectangles_x,rectangles_y,colors); % ACTION!!