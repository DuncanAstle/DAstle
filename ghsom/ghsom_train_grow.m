% GHSOM_TRAIN_GROW  Train one growing SOM, called by ghsom_train.  
%
%
%   sMap = ghsom_train_grow(sData, qe_target, ghTrain, [neighbors, no_edge])
%
%  Input and output arguments ([]'s are optional): 
%   sData     (struct) SOM Toolbox data struct (Training Data)
%   qe_target (scaler) quantisation error target: determines when to stop growing (1..0)
%   ghTrain   (struct) training parameters used
%   neighbors (matrix) the 8 neighbors plus the unit, used for orientation (9 x dims)
%   no_edge   (struct) contains information if there is an edge around (for orientation)
%   sMap      (struct) SOM Toolbox map struct
%   
% For more help, try 'type ghsom_train_grow'.
% See also GHSOM_TRAIN 

%%%%%%%%%%%%% DETAILED DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ghsom_train_grow
% 
% PURPOSE
%
% Grows a Self-Organizing Map used by ghsom_train.
%
% SYNTAX
%
%  sMap = ghsom_train_grow(sD, qe_target, ghTrain, [neighbors, no_edge])
%
% DESCRIPTION
%
% Grows a new GHSOM structure with the given training data
% (SOM Toolbox structure). If the last two arguments are
% omitted it is assumed that the first layer is being trained. 
% The last two arguments are needed to give sublayers an
% orientation. For details on the ghTrain structure passed
% on see ghsom_train.
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
%  sD        (struct) Training data given as a SOM Toolbox data struct.
%                     Use som_data_struct to create it (or have it done by ghsom_train).
%  qe_target (scaler) quantisation error target: determines when to stop growing (1..0)
%                     The value is calculated by ghsom_train using the breadth parameter.
%                     New rows and columns are added until the cummulative or mean (see ghTrain)
%                     quantisation error falls below this value.
%  ghTrain   (struct) training parameters used. 
%                     If no orientation parameters are used then the layer_1st_ ... parameters
%                     are used for training. Otherwise according to the initialization specified
%                     in sub_layer_init will be used. And the first training will be done using
%                     the sub_layer_... parameters. If the quantisation error is yet to high training
%                     will continue (after inserting a row or columne) using the grow_map_...
%                     parameters.
%                     The growth doesnt only depend on the qe_target but also on mbreadth
%                     which limits the number of rows and columns that may be inserted.
%                     Tracking defines how high the verbosity is.
%                     The algorithm setting defines if batch or sequential training is used.
%                     Notice that if the neigh field is set then it will override all the settings
%                     for the different layers.
%  
% OPTIONAL INPUT ARGUMENTS 
%
%  neighbors (matrix) the 8 neighbors plus the unit, used for orientation (9 x dims). 
%                     The whole GHSOM is only defined for a rectangular grid. The units
%                     in the neighborhood are orderd as follows (where 5 is the center):
%                                               1 4 7
%                                               2 5 8
%                                               3 6 9
%                     Non existing neighbors will be zero vectors. To identify them no_edge
%                     is used.
%  no_edge   (struct) contains information if there is an edge around (for orientation)
%                     for example no_edge.right says that all neighbors on the right can
%                     be used because they exist. Notice that different orientation algorithms
%                     handle missing neigbhors differently. It is always assumed that the smallest
%                     map is a 2x2 (it is not possible that no_edge.right == no_edge.left == 0)
%
% OUTPUT ARGUMENTS
% 
%  sMap   (struct) the grown Map in a SOM Toolbox format. 
%
% EXAMPLES
%
% This function has not been intended to call directly. If you would like to create a non-
% hierarchical but growing map, try setting the depth parameter (tau 2) to a value at least
% as high as the width (tau 1) parameter. Specifying the max depth to 1 would do the
% same. If you want to call it directly you would have to create a structure ghTrain and
% pass it as a parameter. See ghsom_train for details.
%
% USED FUNCTIONS FROM SOM TOOLBOX
% 
%  som_bmus, som_randinit, som_lininit, som_batchtrain, som_seqtrain
%
% SEE ALSO
% 
%  ghsom_train  Creates a GHSOM.

% Copyright (c) 2001-2002 by Elias Pampalk, Alvin Chan

% Version 1.0 Alvin Chan 30042001
% Version 2.0 Elias Pampalk 17072001
% Version 2.1 Elias Pampalk 20052002 (fixed some bugs)

function sMap = ghsom_train_grow(sData, qe_target, ghTrain, parent_codebooks, no_edge)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check arguments & initialize

if nargin<3,
    error('Error (ghsom_train): The first three arguments must be specified.');
end

if ghTrain.tracking > 1,
    fprintf('   Growing SOM: ');
    if nargin==3,
        fprintf('Toplayer ');
    else
        fprintf('Sublayer ');
    end
    fprintf('-- Target QE: %3.2f\n',qe_target);
end

% choose settings for the first training
if nargin == 3,
    r0     = ghTrain.layer_1st_radius;
    neigh0 = ghTrain.layer_1st_neigh;
    init0  = ghTrain.layer_1st_init;
    alpha0 = 0.5; % only needed if seq, rough training
else
    r0     = ghTrain.sub_layer_radius;
    neigh0 = ghTrain.sub_layer_neigh;
    init0  = ghTrain.sub_layer_init;
    alpha0 = 0.05; % only needed if seq, finetunning
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% action

switch init0,
 case 'randinit', sMap = som_randinit(sData,'msize',[2 2],'rect','sheet');
 case 'lininit',  sMap = som_lininit(sData,'msize',[2 2],'rect','sheet');
 case 'sample', sMap = ghsom_train_grow_sampleinit(sData, 'msize',[2 2],'rect','sheet');
 otherwise, sMap = ghsom_train_grow_orientation(parent_codebooks, no_edge, init0);
end

if ghTrain.tracking > 1,
    [bmus,qerr] = som_bmus(sMap,sData);
    sqerr=sum(qerr);
    fprintf('       After initialization sum(qerr): %f\n',sqerr);
end    

switch ghTrain.algorithm,
 case 'batch',
    sMap = som_batchtrain(sMap,sData,'tracking',0,'radius',r0,'neigh',neigh0);
 case 'seq',
    sMap = som_seqtrain(sMap,sData,'tracking',0,'radius',r0,'alpha_ini',alpha0, ...
        'alpha_type','inv','neigh',neigh0);
 otherwise, error('Error (ghsom_train_grow): unknown algorithm type!');
end

if ghTrain.tracking > 1,
    [bmus,qerr] = som_bmus(sMap,sData);
    sqerr2=sum(qerr);
    fprintf('       After training (2x2 map) sum(qerr): %f\n',sqerr2);
    if sqerr == sqerr2,
        disp('      WARNING!!! No improvement through training!!\n\n');
    end
end    

[dummy qerrs] = som_bmus(sMap,sData);
switch ghTrain.qe_type, 
 case 'cummulative', qe = sum(qerrs);
 case 'mean',        qe = mean(qerrs);
 otherwise, error('Error (ghsom_train_grow): unknown quantisation error type!');
end

insertions=0;
% continue growing until the quantisation error is low enough
while qe > qe_target & insertions<ghTrain.mbreadth,  
    insertions=insertions+1;
    [unit_mqe, unit_mdist_neigbhor] = find_mconflict_units(sMap,sData,ghTrain.qe_type);        
    sMap = add_row_or_col(sMap,unit_mqe,unit_mdist_neigbhor,ghTrain.tracking);
                
    switch ghTrain.algorithm,
     case 'batch', sMap = som_batchtrain( ...
             sMap,sData,'tracking',0,'radius',ghTrain.grow_map_radius,...
             'neigh',ghTrain.grow_map_neigh);
     case 'seq',   sMap = som_seqtrain( ...
             sMap,sData,'tracking',0,'radius',ghTrain.grow_map_radius,'alpha_ini',0.05,...
             'alpha_type','inv', 'neigh',ghTrain.grow_map_neigh);
    end
        
    [dummy qerrs] = som_bmus(sMap,sData);
    switch ghTrain.qe_type, 
     case 'cummulative', qe = sum(qerrs);
     case 'mean',        qe = mean(qerrs);
     otherwise, error('Error (ghsom_train_grow): unknown quantisation error type!');
    end
       
    if ghTrain.tracking > 1,
        fprintf(' - QE after Insertion: %3.2f\n',qe);
    end
end
if insertions==ghTrain.mbreadth,
    disp('Warning: Max breadth reached!');
end
%% END Function: ghsom_train_grow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sub Functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [unit_mqe, unit_mdist_neigbhor] = ...
    find_mconflict_units(sMap,sData,qe_type)
%
% Calculates Cummulative qerr for each unit compared with input vectors
% Returns unit with the highest quantisation error.
%
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

[qe_max, unit_mqe] = max(unit_qerr);

% NOTE: som_eucdist2(sMap,sMap) would ignore mask

% find neighbor to unit with mqe which is most furthest away
rows = sMap.topol.msize(1);
cols = sMap.topol.msize(2);

mdist = -1; % maximum distance to immediate neighbors
unit_mdist_neigbhor  = -1; % index of this neigbhor
if unit_mqe-rows > 0, % left unit exists
    unit_mdist_neigbhor = unit_mqe - rows;
    mdist = sqrt(sum((sMap.mask'.*(sMap.codebook(unit_mqe,:)-sMap.codebook(unit_mqe-rows,:))).^2));
    %mdist = sqrt(sum((sMap.mask.*(sMap.codebook(unit_mqe)-sMap.codebook(unit_mqe-rows))).^2));
end
if unit_mqe+rows <= rows*cols, % right unit exists
    dist = sqrt(sum((sMap.mask'.*(sMap.codebook(unit_mqe,:)-sMap.codebook(unit_mqe+rows,:))).^2));
    %dist = sqrt(sum((sMap.mask.*(sMap.codebook(unit_mqe)-sMap.codebook(unit_mqe+rows))).^2));
    if dist > mdist,
        mdist = dist;
        unit_mdist_neigbhor = unit_mqe+rows;
    end
end
if mod(unit_mqe,rows) ~= 1, % upper unit exists
    dist = sqrt(sum((sMap.mask'.*(sMap.codebook(unit_mqe,:)-sMap.codebook(unit_mqe-1,:))).^2));
    %dist = sqrt(sum((sMap.mask.*(sMap.codebook(unit_mqe)-sMap.codebook(unit_mqe-1))).^2));
    if dist > mdist,
        mdist = dist;
        unit_mdist_neigbhor = unit_mqe-1;
    end
end
if mod(unit_mqe,rows) ~= 0, % lower unit exists
    dist = sqrt(sum((sMap.mask'.*(sMap.codebook(unit_mqe,:)-sMap.codebook(unit_mqe+1,:))).^2));
    %dist = sqrt(sum((sMap.mask.*(sMap.codebook(unit_mqe)-sMap.codebook(unit_mqe+1))).^2));
    if dist > mdist,
        mdist = dist;
        unit_mdist_neigbhor = unit_mqe+1;
    end
end
%% END Function: find_mconflict_units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sMap = add_row_or_col(sMap,unit1,unit2,tracking)
%
% Returns sMap after growing one step. this is done by adding a row
% (if the units are on the same columne) or adding a new column
% (if the units are on the same row) between these units.

rows = sMap.topol.msize(1);
cols = sMap.topol.msize(2);

TRANSPOSE=0;
if abs(unit1-unit2)==1 % add a new row
    TRANSPOSE=1;
    idx = rot90(reshape(1:rows*cols,rows,cols));
    [dummy i_idx] = sort(idx(:));
    sMap.codebook = sMap.codebook(idx(:),:); 
    unit1 = i_idx(unit1);
    unit2 = i_idx(unit2);
    if tracking>1, fprintf('           Adding Row (%d,%d)',rows+1,cols); end    
    sMap.topol.msize(1) = rows+1;
    cols = rows;
    rows = sMap.topol.msize(2);
else % add a new column
    if tracking>1, fprintf('           Adding Column (%d,%d)',rows,cols+1); end
    sMap.topol.msize(2) = cols+1;    
end

c1 = ceil(unit1/rows);
c2 = ceil(unit2/rows);

mr=min([c1,c2])*rows;
new_codebook = [...
        sMap.codebook(1:mr,:); ...
        1/2*(sMap.codebook(((1:rows)+mr)-rows,:)+sMap.codebook((1:rows)+mr,:)); ...
        sMap.codebook((mr+1):end,:)];

if TRANSPOSE, % row, added, transpose back
    idx = rot90(rot90(rot90(reshape(1:rows*(cols+1),rows,cols+1))));
    sMap.codebook = new_codebook(idx(:),:);
else % col, added, dont need to transpose back
    sMap.codebook = new_codebook;
end
sMap.labels = cell(size(sMap.codebook,1),1); % Add extra label fields for new units
%% END Function: add_row_or_col
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
