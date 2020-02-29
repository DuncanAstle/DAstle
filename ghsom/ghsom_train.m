% GHSOM_TRAIN  Train (create) a Growing Hierarchical SOM 
%
%  ghMap = ghsom_train(D, [argID, value, ...])
%
%  ghMap = ghsom_train(D)
%  [ghMap ghTrain] = ghsom_train(D, 'breadth', 0.3, 'depth', 0.01)
%
%  Input and output arguments ([]'s are optional): 
%   D        (matrix) training data, size dlen x dim
%            (struct) SOM Toolbox data struct
%   [argID,  (string) See below.
%    value]  (varies) 
%
%   ghMap    (struct) ghrowing hierarchical map struct
%   ghTrain  (struct) training parameters used
%
% Here are the valid argument IDs and corresponding values.
%   'breadth'               (scaler) Range from 1 to 0, controlls breadth of maps 
%   'depth'                 (scaler) Range from 1 to 0, controlls depth of GHSOM
%   'depth_dataitems'       (scaler) Units with less dataitems cannot be expanded
%   'layer_1st_radius'      (vector) neighborhood radiuses, length 1, 2 or trainlen
%   'layer_1st_neigh'       (string) neighborhood function: 'gaussian', 'cutgauss',
%                                    'ep', or 'bubble'
%   'layer_1st_init'        (string) initialization: 'randinit' (default), 'sample', or 'lininit' 
%   'layer_1st_trainlen'    (scaler) training lengths in epochs
%   'sub_layer_radius'      (vector) neighborhood radiuses, length 1, 2 or trainlen
%   'sub_layer_neigh'       (string) neighborhood function: 'gaussian', 'cutgauss',
%                                    'ep', or 'bubble'
%   'sub_layer_init'        (string) orientation: 'dittenbach', 'mirror' (default), 'sample', 
%                                    'randinit', or 'lininit'
%   'sub_layer_trainlen'    (scaler) training lengths in epochs
%   'grow_map_radius'       (vector) neighborhood radiuses, length 1, 2 or trainlen
%   'grow_map_radius_msize' (string) radius relation to mapsize: 'none' or 'linear' 
%   'grow_map_neigh'        (string) neighborhood function: 'gaussian', 'cutgauss',
%                                    'ep', or 'bubble'
%   'grow_map_trainlen'     (scaler) training lengths in epochs
%   'neigh'                 (string) neighborhood function, overwrites others
%   'algorithm'             (string) which algorithm to use: 'batch' (default) or 'seq'
%   'qe_type'               (string) quantisation error: 'mean' or 'cummulative' (default) 
%   'mbreadth'              (scaler) escape infinite loops in width
%   'mdepth'                (scaler) escape infinite loops in depth 
%   'tracking'              (scaler) how much to report, levels 0-3, default == 1
%   
% For more help, try 'type ghsom_train'.
% See also SOM_DATA_STRUCT, GHSOM_TRAIN_GROW 

%%%%%%%%%%%%% DETAILED DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ghsom_train
% 
% PURPOSE
%
% Creates a trained Growing Hierarchical Self-Organizing Map.
%
% SYNTAX
%
%  ghMap = ghsom_train(D)
%  ghMap = ghsom_train(sD)
%  [ghMap ghTrain] = ghsom_train(...,'argID',value,...);
%
% DESCRIPTION
%
% Creates a new GHSOM structure with the given training data
% (sD or D).  If no optional arguments (argID, value) are given, 
% a default training is done. Using optional arguments the 
% training parameters can be specified. Returns the newly created 
% and trained GHSOM and the parameters used.
%
% REFERENCES
%
% Dittenbach, M., Merkl, D., and Rauber, A. "The Growing Hierarchical Self-Organizing Map"
%    S. Amari and C. L. Giles and M. Gori and V. Puri, editors, Proceedings of the 
%    International Joint Conference on Neural Networks (IJCNN 2000), vol. 6, pp. 15-19, 
%    July 24-27, 2000, Como, Italy, IEEE Computer Society.
%
% REQUIRED INPUT ARGUMENTS
%
%  D           Training data.
%     (struct) data struct
%     (matrix) data matrix, size [dlen dim]
%  
% OPTIONAL INPUT ARGUMENTS 
%
%  argID (string) Argument identifier string (see below).
%  value (varies) Value for the argument (see below).
%
%  The optional arguments can be given as 'argID',value -pairs. If an
%  argument is given value multiple times, the last one is
%  used. The valid IDs and corresponding values are listed below.
%
%  Below is the list of valid arguments: 
%   'breadth'               (scaler) Range from 1 to 0, controlls breadth of maps
%                                    default is 0.3, published in papers as tau_1
%                                    setting breadth to 1 would lead to only 2x2 maps
%                                    setting it to 0 would lead to a (too) huge map
%                                    the map will grown in width and height until
%                                    breadth * cummulative QE < cummulative QE of the
%                                    upper map (QE = quantisation error)
%                                    (the upper layer of the 1st layer is a 
%                                    single unit map)
%   'depth'                 (scaler) Range from 1 to 0, controlls depth of GHSOM
%                                    default is 0.03, published in papers as tau_2
%                                    a new map is created for the data mapped onto one
%                                    unit if the quantisation error (QE) of the unit
%                                    is bigger than depth * cummulative QE.
%                                    only values smaller than breadth make sense since
%                                    the QE of a unit will always be smaller as the QE
%                                    of all units on the map. 
%                                    setting depth to 1 will lead to no hierarchy,
%                                    setting depth to 0 will lead to a very deep branches.
%   'depth_dataitems'       (scaler) Units with less dataitems cannot be expanded. Default is
%                                    set to 4, so if a unit has only (or less then) 4 dataitems 
%                                    mapped regardless is quantisation error, it will not 
%                                    be furter expanded.
%   'layer_1st_radius'      (vector) neighborhood radiuses, length 1, 2 or trainlen
%                                    length = 1: radius = linspace(radius,radius/4,trainlen)
%                                    length = 2: radius = linspace(radius(1),radius(2),trainlen)
%                                    length > 2: the vector given neighborhood
%                                       radius for each step separately
%                                       trainlen = length(radius)
%   'layer_1st_neigh'       (string) neighborhood function: 'gaussian' (default), 'cutgauss',
%                                    'ep' or 'bubble'
%   'layer_1st_init'        (string) initialization: 'randinit' (default), 'sample', or 'lininit' 
%   'layer_1st_trainlen'    (scaler) training lengths in epochs
%   'sub_layer_radius'      (vector) neighborhood radiuses, length 1, 2 or trainlen
%                                    length = 1: radius = linspace(radius,radius/4,trainlen)
%                                    length = 2: radius = linspace(radius(1),radius(2),trainlen)
%                                    length > 2: the vector given neighborhood
%                                       radius for each step separately
%                                       trainlen = length(radius)
%   'sub_layer_neigh'       (string) neighborhood function: 'gaussian' (default), 'cutgauss',
%                                    'ep' or 'bubble'
%   'sub_layer_init'        (string) orientation: 'dittenbach', 'mirror' (default), 'sample'
%                                    (which randomly selects dataitems), 'randinit', or 'lininit'
%   'sub_layer_trainlen'    (scaler) training lengths in epochs
%   'grow_map_radius'       (vector) neighborhood radiuses, length 1, 2 or trainlen
%                                    length = 1: radius = linspace(radius,radius/4,trainlen)
%                                    length = 2: radius = linspace(radius(1),radius(2),trainlen)
%                                    length > 2: the vector given neighborhood
%                                       radius for each step separately
%                                       trainlen = length(radius)
%   'grow_map_radius_msize' (string) radius relation to mapsize: 'none' or 'linear' 
%   'grow_map_neigh'        (string) neighborhood function: 'gaussian' (default), 'cutgauss',
%                                    'ep' or 'bubble'
%   'grow_map_trainlen'     (scaler) training lengths in epochs
%   'neigh'                 (string) neighborhood function, overwrites the individual settings
%   'algorithm'             (string) which algorithm to use: 'batch' (default) or 'seq'
%   'qe_type'               (string) quantisation error: 'mean' or 'cummulative' (default)
%                                    this referes to the tau_1 and tau_2 limits. shall they be
%                                    calculated based on the mean or rather on the cummalative
%                                    quantisation error...
%   'mbreadth'           (scaler) escape infinite loops in width, it might occur that the
%                                    algorithm keeps adding rows and columns to a map. To avoid this
%                                    this value is set to 20 by default. A warning occurs if this
%                                    happens.
%   'mdepth'             (scaler) escape infinite loops in depth, it might occur that the
%                                    GHSOM keeps growing in depth. To avoid this this value
%                                    is set to 100 by default. A warning occurs if this happens.
%   'tracking'   (scalar) tracking level: silent 0, 1 (default), 2 or 3 maximum verbosity
%
% OUTPUT ARGUMENTS
% 
%  ghMap   (struct) the created and trained growing hierarchical map
%  ghTrain (struct) the parameters used for training the GHSOM
%
% EXAMPLES
%
% Simplest case:
%  ghM = ghsom_train(D);  
%  ghM = ghsom_train(sD);  
%
% To change the tracking level, 'tracking' argument is specified:
%  ghM = ghsom_train(D,'tracking',3);
%
% To change training parameters, the optional arguments 'breadth', 'depth', 
% 'depth_dataitems', 'layer_1st_radius', 'layer_1st_neigh', 'layer_1st_init', 
% 'layer_1st_trainlen', 'sub_layer_radius', 'sub_layer_neigh', 'sub_layer_init', 
% 'sub_layer_trainlen', 'grow_map_radius', 'grow_map_radius_msize', 'grow_map_neigh', 
% 'grow_map_trainlen', 'neigh', 'algorithm', 'qe_type', 'mbreadth', and 
% 'mdepth' are used.
% 'layer_1st' is the first 2x2 Map which is trained. 'sub_layer' refers to any new
% created 2x2 maps other than the 1st layer. 'grow_map' parameters are used after a
% new columne or row has been added to a map.
%  ghMap = ghsom_train(D,'breadth',0.3,'depth',0.03,'orient','none',...
%  'neigh','gaussian','layer_1st_init','lininit')
%
% Notice that if the sequential algorithm is used than the SOM Toolbox default value
% for rough training is used on the first map (alpha_ini = 0.5, alpha_type = 'inv')
% and for the 'sub_layer' and the 'grow_map' the finetuning values are used 
% (alpha_ini = 0.05, alpha_type = 'inv').
%
% (DIRECTLY) USED FUNCTIONS FROM SOM TOOLBOX
% 
%  som_divide, som_bmus, som_data_struct, som_vis_coords
%
% SEE ALSO
% 
%  som_data_struct   Create a SOM data struct (SOM Toolbox).
%  ghsom_train_grow  Grows (adding rows and columns) one SOM.

% Copyright (c) 2001 by Elias Pampalk, Alvin Chan

% Version 1.0 Alvin Chan 30042001
% Version 2.0 Elias Pampalk 17072001

function [ghMap, ghTrain] = ghsom_train(sData, varargin) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check arguments & initialize

new_growing = 0;

if nargin<1,
    error('Error (ghsom_train): Input argument sData not defined.');
end

% sData
if ~isstruct(sData),
    sData = som_data_struct(sData,'name',inputname(1)); 
    
end
[dlen dim] = size(sData.data);

% varargin
% parameters default set
ghTrain.breadth=0.3;
ghTrain.depth=0.03;
ghTrain.depth_dataitems=4;
ghTrain.layer_1st_radius=[1.5,0.1];
%ghTrain.layer_1st_radius=[1.5,0.5];
ghTrain.layer_1st_neigh='gaussian';
ghTrain.layer_1st_init='lininit';
%ghTrain.layer_1st_init='randinit';
ghTrain.layer_1st_trainlen=20;
ghTrain.sub_layer_radius=[0.5,0.1];
ghTrain.sub_layer_neigh='gaussian';
ghTrain.sub_layer_init='mirror';
%ghTrain.sub_layer_trainlen=10;
ghTrain.sub_layer_trainlen=35;
ghTrain.grow_map_radius=[0.8,0.1];
ghTrain.grow_map_radius_msize='none';
ghTrain.grow_map_neigh='gaussian';
ghTrain.grow_map_trainlen=15;
ghTrain.neigh=NaN;
ghTrain.algorithm='batch';
ghTrain.qe_type='cummulative';
ghTrain.mbreadth=20;
ghTrain.mdepth=100;
ghTrain.tracking=1;

i=1; 
while i<=length(varargin), 
    argok = 1; 
    if ischar(varargin{i}), 
        switch varargin{i}, 
            % argument IDs
        case 'breadth', i=i+1; ghTrain.breadth = varargin{i}; 
        case 'depth', i=i+1; ghTrain.depth = varargin{i};
        case 'depth_dataitems', i=i+1; ghTrain.depth_dataitems = varargin{i};
        case 'layer_1st_radius', 
            i=i+1; 
            l = length(varargin{i}); 
            if l==1, 
                ghTrain.layer_1st_radius = [varargin{i},varargin{i}/4]; 
            elseif l==2 
                ghTrain.layer_1st_radius = [varargin{i}(1),varargin{i}(2)];  
            elseif l>2
                ghTrain.layer_1st_radius = varargin{i};
            end
        case 'layer_1st_neigh', i=i+1; ghTrain.layer_1st_neigh = varargin{i};
        case 'layer_1st_init', i=i+1; ghTrain.layer_1st_init = varargin{i};
        case 'layer_1st_trainlen', i=i+1; ghTrain.layer_1st_trainlen = varargin{i};
        case 'sub_layer_radius', i=i+1; ghTrain.sub_layer_radius = varargin{i};
            i=i+1; 
            l = length(varargin{i}); 
            if l==1, 
                ghTrain.sub_layer_radius = [varargin{i},varargin{i}/4]; 
            elseif l==2 
                ghTrain.sub_layer_radius = [varargin{i}(1),varargin{i}(2)];  
            elseif l>2
                ghTrain.sub_layer_radius = varargin{i};
            end
        case 'sub_layer_neigh', i=i+1; ghTrain.sub_layer_neigh = varargin{i};
        case 'sub_layer_init', i=i+1; ghTrain.sub_layer_init = varargin{i};
        case 'sub_layer_trainlen', i=i+1; ghTrain.sub_layer_trainlen = varargin{i};
        case 'grow_map_radius', 
            i=i+1; 
            l = length(varargin{i}); 
            if l==1, 
                ghTrain.grow_map_radius = [varargin{i},varargin{i}/4]; 
            elseif l==2 
                ghTrain.grow_map_radius = [varargin{i}(1),varargin{i}(2)];  
            elseif l>2
                ghTrain.grow_map_radius = varargin{i};
            end
        case 'grow_map_radius_msize', i=i+1; ghTrain.grow_map_radius_msize = varargin{i};
        case 'grow_map_neigh', i=i+1; ghTrain.grow_map_neigh = varargin{i};
        case 'grow_map_trainlen', i=i+1; ghTrain.grow_map_trainlen = varargin{i};
        case 'neigh', i=i+1; ghTrain.neigh = varargin{i};
        case 'algorithm', i=i+1; ghTrain.algorithm = varargin{i};
        case 'qe_type', i=i+1; ghTrain.qe_type = varargin{i};
        case 'mbreadth', i=i+1; ghTrain.mbreadth = varargin{i};
        case 'mdepth', i=i+1; ghTrain.mdepth = varargin{i};
        case 'tracking', i=i+1; ghTrain.tracking = varargin{i};
        otherwise argok=0; 
        end
    else argok=0;
    end
    if ~argok, 
        disp(['WARNING: (ghsom_train) Ignoring invalid argument #' num2str(i+2)]); 
    end
    i = i+1; 
end


% radius
if length(ghTrain.layer_1st_radius)>2
    ghTrain.layer_1st_trainlen=lenght(ghTrain.layer_1st_radius);
else
    if strcmp(ghTrain.algorithm,'batch'),
        ghTrain.layer_1st_radius=linspace(ghTrain.layer_1st_radius(1),...
            ghTrain.layer_1st_radius(2),ghTrain.layer_1st_trainlen);
    else % if seq then multiply by dataitems (dont use epochs)
        ghTrain.layer_1st_radius=linspace(ghTrain.layer_1st_radius(1),...
            ghTrain.layer_1st_radius(2),ghTrain.layer_1st_trainlen * dlen);
    end
end
if length(ghTrain.sub_layer_radius)>2
    ghTrain.sub_layer_trainlen=lenght(ghTrain.sub_layer_radius);
else
    if strcmp(ghTrain.algorithm,'batch'),
        ghTrain.sub_layer_radius=linspace(ghTrain.sub_layer_radius(1),...
            ghTrain.sub_layer_radius(2),ghTrain.sub_layer_trainlen);
    else % if seq then multiply by dataitems (dont use epochs)
        ghTrain.sub_layer_radius=linspace(ghTrain.sub_layer_radius(1),...
            ghTrain.sub_layer_radius(2),ghTrain.sub_layer_trainlen * dlen);
    end
end
if length(ghTrain.grow_map_radius)>2
    ghTrain.grow_map_trainlen=lenght(ghTrain.grow_map_radius);
else
    if strcmp(ghTrain.algorithm,'batch'),
        ghTrain.grow_map_radius=linspace(ghTrain.grow_map_radius(1),...
            ghTrain.grow_map_radius(2),ghTrain.grow_map_trainlen);
    else % if seq then multiply by dataitems (dont use epochs)
        ghTrain.grow_map_radius=linspace(ghTrain.grow_map_radius(1),...
            ghTrain.grow_map_radius(2),ghTrain.grow_map_trainlen * dlen);
    end
end

% override neighborhood settings if neigh is defined
if ~isnan(ghTrain.neigh),
    ghTrain.layer_1st_neigh = ghTrain.neigh;
    ghTrain.sub_layer_neigh = ghTrain.neigh;
    ghTrain.grow_map_neigh = ghTrain.neigh;
end

if strcmp(ghTrain.layer_1st_init,'mirror') | strcmp(ghTrain.layer_1st_init,'dittenbach'),
    error('Error (ghsom_train): Cannot use orientation initialization on first layer!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% action

% ghMap structure:
%   layer             (vector) layer on which map i is.
%   parent            (vector) parent map from which map i descended.
%   parent_unit       (vector) parent map unit from which map i descended.
%   sMap                (cell) the sMap (SOM Toolbox) of the (sub)map i.
%   expand_units        (cell) which units are expanded on map i.
%   expand_units_qe     (cell) qerr of the these units.
%   dataitems           (cell) indexes of dataitems belonging to map i.

% (theoretical) layer 0 - one unit (neuron) map
[bmus,qerrs] = som_bmus(mean(sData.data),sData);

switch ghTrain.qe_type, 
 case 'cummulative', qe0 = sum(qerrs(:,1));
 case 'mean',        qe0 = mean(qerrs(:,1));
 otherwise, error('Error (ghsom_train): unknown quantisation error type!');
end
qe_depth   = qe0 * ghTrain.depth; % only calculated once

if ghTrain.tracking>0,
    my_time=clock;
    my_cpu_time0=cputime;
    fprintf('Starting GHSOM Training at %2d:%2d:%2d (hh:mm:ss)\n',my_time(4),my_time(5),round(my_time(6)));
    if ghTrain.tracking>1
        fprintf('Layer 0\n');
        fprintf('  QE: %3.2f\n', qe0);
        fprintf('  Breadth QE: %3.2f\n', qe0 * ghTrain.breadth);
        fprintf('  Depth QE: %3.2f\n\n', qe_depth);
    end
end

% special case: 1st layer
ghMap.layer(1) = 1;
if new_growing,
    ghMap.sMap{1} = ghsom_train_grow2(sData);
else
    ghMap.sMap{1} = ghsom_train_grow(sData, qe0 * ghTrain.breadth, ghTrain); % ***
end
ghMap.dataitems{1}=1:dlen;
ghMap.parent(1) = 0;
ghMap.parent_unit(1) = 0;

parentIdx = 0; 

% main loop - loop through all parents, 
% end when all parents have been seen 
% parents can have children, which are
% added to the end of the 'todo' list
while parentIdx < length(ghMap.sMap), 
    parentIdx = parentIdx + 1;
            
    if ghMap.layer(parentIdx)<ghTrain.mdepth,
        % on parent map find units to expand and their qerr
        [ghMap.expand_units{parentIdx}, ghMap.expand_units_qe{parentIdx}] = find_expand_units( ...
            ghMap.sMap{parentIdx}, ...
            ghsom_idx2data(sData,ghMap.dataitems{parentIdx}), ...
            qe_depth, ...
            ghTrain.depth_dataitems, ...
            ghTrain.qe_type, ...
            ghTrain.tracking);
    else % too far down already
        disp('Warning: Max depth reached!');
        ghMap.expand_units{parentIdx} = [];
        ghMap.expand_units_qe{parentIdx} = [];
    end
    
    if ghTrain.tracking > 1,
        [dummy qerrs] = som_bmus(ghMap.sMap{parentIdx},ghsom_idx2data(sData,ghMap.dataitems{parentIdx}));
        switch ghTrain.qe_type, 
         case 'cummulative', qe = sum(qerrs);
         case 'mean',        qe = mean(qerrs);
        end
        fprintf('   Layer: %d, Parent: %2d, Parent Unit: %2d, QE: %3.2f, Children: %d\n', ...
            ghMap.layer(parentIdx), ...
            ghMap.parent(parentIdx), ...
            ghMap.parent_unit(parentIdx), ...
            qe, ...
            length(ghMap.expand_units{parentIdx}));
    end
    
    % loop through all children of one parent and create their maps
    for i=1:length(ghMap.expand_units{parentIdx}), % if there is none, expand_units{.} = []
        thisMapIdx = length(ghMap.sMap) + 1; % add this child to end of list
        ghMap.layer(thisMapIdx) = ghMap.layer(parentIdx) + 1; % immediate below parent
        
        % find data belonging to this child
        [dummy ghMap.dataitems{thisMapIdx}] = som_divide( ...
            ghMap.sMap{parentIdx}, ...
            ghsom_idx2data(sData,ghMap.dataitems{parentIdx}), ... 
            ghMap.expand_units{parentIdx}(i), ...
            'index');
        ghMap.dataitems{thisMapIdx} = ghMap.dataitems{parentIdx}(ghMap.dataitems{thisMapIdx});
        
        % get the 8 adjectant neigbhors plus parent unit itself 
        [neighbors no_edge] = get_neighbors( ...
            ghMap.sMap{parentIdx}, ...
            ghMap.expand_units{parentIdx}(i), ...
            ghTrain.tracking);
        
        % create child map
        if new_growing,
            ghMap.sMap{thisMapIdx} = ghsom_train_grow2(ghsom_idx2data(sData, ghMap.dataitems{thisMapIdx}),neighbors,no_edge);
        else
            ghMap.sMap{thisMapIdx} = ghsom_train_grow( ... % ***
                ghsom_idx2data(sData, ghMap.dataitems{thisMapIdx}), ...
                ghTrain.breadth * ghMap.expand_units_qe{parentIdx}(i), ...
                ghTrain, ...    
                neighbors, ...
                no_edge);
        end
        
        ghMap.parent(thisMapIdx) = parentIdx;
        ghMap.parent_unit(thisMapIdx) = ghMap.expand_units{parentIdx}(i);
        
        if ghTrain.tracking > 2,
            fprintf('  Created Child:\n')
            fprintf('    Layer: %d, Parent: %2d, Parent Unit: %2d\n', ...
                ghMap.layer(thisMapIdx), ...
                ghMap.parent(thisMapIdx), ...
                ghMap.parent_unit(thisMapIdx));
        end
    end
end % main while loop

if ghTrain.tracking>0,
    my_time=clock;
    my_cpu_time_total=cputime-my_cpu_time0;
    fprintf('Ended GHSOM Training at %2d:%2d:%2d (hh:mm:ss)\n',my_time(4),my_time(5),round(my_time(6)));
    fprintf('CPU time consumed: %dm %ds\n', floor(my_cpu_time_total/60), round(my_cpu_time_total-60*floor(my_cpu_time_total/60)));
end

%% END Function: ghsom_train
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sub Functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [neighbors, no_edge] = get_neighbors(sMap, unit, tracking)
%
% find 8 neighbors of a unit (rectangular grid)
% if unit is not available all zeros in the dims
% in struct no_edge it says weather there is 
% no edge (up, down, left or right). if there is
% no edge it means that there are neighor units in
% that area.
% neighbors is stored in SOM Toolbox numbering scheme
%   1 4 7
%   2 5 8
%   3 6 9
% where 5 ist the center unit, always the unit which
% is being expanded. so if there is an upper edge
% and there are no units above the expanded unit,
% then the neighbors 2,3,5,6,8,9 will be returned.

rows = sMap.topol.msize(1);
cols = sMap.topol.msize(2);

coords = som_vis_coords('rect',sMap.topol.msize);
unit_coords = coords(unit,:);

neighbors = zeros(9,size(sMap.codebook,2));

no_edge.up    = unit_coords(2)-1>=1;
no_edge.down  = unit_coords(2)+1<=rows;
no_edge.left  = unit_coords(1)-1>=1;
no_edge.right = unit_coords(1)+1<=cols;

if no_edge.left & no_edge.up,
    neighbors(1,:) = sMap.codebook(unit-rows-1,:);  % up-left
end        
if no_edge.left, 
    neighbors(2,:) = sMap.codebook(unit-rows,:);    % center-left
end        
if no_edge.left & no_edge.down, 
    neighbors(3,:) = sMap.codebook(unit-rows+1,:);  % down-left
end        
if no_edge.up, 
    neighbors(4,:) = sMap.codebook(unit-1,:);       % up-center
end               
neighbors(5,:) = sMap.codebook(unit,:);             % center (always exists)  
if no_edge.down, 
    neighbors(6,:) = sMap.codebook(unit+1,:);       % down-center
end           
if no_edge.up & no_edge.right, 
    neighbors(7,:) = sMap.codebook(unit+rows-1,:);  % up-right
end 
if no_edge.right, 
    neighbors(8,:) = sMap.codebook(unit+rows,:);    % center-right
end 
if no_edge.right & no_edge.down, 
    neighbors(9,:) = sMap.codebook(unit+rows+1,:);  % down-right
end

if tracking==3,
    fprintf('    Got neighbors for unit (%d,%d).\n',unit_coords(1),unit_coords(2));
    if ~no_edge.up,
        fprintf('        Edge above exists.\n');
    end
    if ~no_edge.down,
        fprintf('        Edge below exists.\n');
    end
    if ~no_edge.left,
        fprintf('        Edge left exists.\n');
    end
    if ~no_edge.right,
        fprintf('        Edge right exists.\n');
    end
end
%% END Function: get_neighbors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [expand_units, expand_units_qe] = find_expand_units( ...
    sMap, sData, qe_depth, depth_dataitems, qe_type, tracking)
%
% find the units from a map which satisfy the criteria for expandation.
% returns the units and their QE
% the criteria depend are:
%   1. qe of unit > 'qe_depth' 
%   2. dataitems on unit > 'depth_dataitem'

[bmus, qerr] = som_bmus(sMap, sData);
munits = size(sMap.codebook,1);
unit_qerr = zeros(munits,1);
unit_mapped = zeros(munits,1);

for i = 1:munits,
    unit_qerr(i) = sum(qerr(find(bmus==i)));
    unit_mapped(i) = sum(bmus==i);
end
switch qe_type, 
 case 'mean',
    idx = find(unit_mapped>0);
    unit_qerr(idx) = unit_qerr(idx) ./ unit_mapped(idx);
end
 
if tracking==3,
    fprintf('   Looking for units to expand ...\n');
    fprintf('       Need at least QE: %f\n', qe_depth);
    fprintf('       Units satisfying "depth" criteria:');
    fprintf(' %d',find(unit_qerr>qe_depth));
    if isempty(find(unit_qerr>qe_depth)),
        fprintf('none');
    end
    fprintf('.\n');
    fprintf('       Units satisfying "depth_dataitems" criteria:');
    fprintf(' %d',find(unit_mapped>depth_dataitems));
    if isempty(find(unit_mapped>depth_dataitems)),
        fprintf('none');
    end
    fprintf('.\n');
end

idx = find(unit_qerr>qe_depth & unit_mapped>depth_dataitems);
expand_units = idx;
expand_units_qe = unit_qerr(idx);
if tracking==3,
    if isempty(idx),
        fprintf('     Found none.\n');
    else
        fprintf('     Found: Unit %2d, QE: %2.2f, Dataitems: %3d.\n',[idx,expand_units_qe,unit_mapped(idx)]');
    end
end
%% END Function: find_expand_units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%