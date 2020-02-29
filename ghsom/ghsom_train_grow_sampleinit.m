% GHSOM_TRAIN_GROW_SAMPLEINIT  Initialize a SOM by using random samples, called by ghsom_train_grow.  
%
%
%   sMap = ghsom_train_grow_sampleinit(sData, msize)
%
%   sData     (struct) SOM Toolbox data struct (Training Data)
%   msize     (vector) map size: [rows cols] 
%   sMap      (struct) SOM Toolbox map struct
%   
% For more help, try 'type ghsom_train_grow_sampleinit'.
% See also GHSOM_TRAIN_GROW

%%%%%%%%%%%%% DETAILED DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ghsom_train_grow_sampleinit
% 
% PURPOSE
%
% Initialize a SOM (usually a 2x2, the first layer of the GHSOM) using randomly selected dataitems.
%
% SYNTAX
%
%  sM = ghsom_train_grow_sampleinit(sD, msize)
%  sM = ghsom_train_grow_sampleinit(sD, [2 2])
%
% DESCRIPTION
%
% Randomly selects elements out of 'sD.data'. If there are more dataitems than map units
% no dataitem will be used more than once. Otherwise if there are more map units than
% dataitems a dataitem can be theoretically used as many times as there are map units.
%
% REQUIRED INPUT ARGUMENTS
%
%  sD        (struct) Training data given as a SOM Toolbox data struct.
%                     Use som_data_struct to create it (or have it done by ghsom_train).
%  msize     (vector) Map size: msize(1) are the number of rows, msize(2) are the number of columns.
%                     Notice that only rectangular sheet type maps are supported.
%
% OUTPUT ARGUMENTS
% 
%  sMap   (struct) the grown Map in a SOM Toolbox format. 
%
% SEE ALSO
% 
%  ghsom_train_grow  Grows (adding rows and columns) one SOM.

% Copyright (c) 2001 by Elias Pampalk

% Version 1.0 Elias Pampalk 17072001

function sM = ghsom_train_grow_sampleinit(sD, msize)

munits = prod(msize);
[dlen dims] = size(sD.data);

if dlen < munits,
    fprintf('Warning: sample initializing SOM with %d units using only %d vectors!',munits,dlen);
end

sMap = som_map_struct(dims, 'msize', msize, 'sheet', 'rect'); 

if dlen > munits,
    idx = randperm(dlen);
else
    idx = ceil(rand(1,3)*5);
end
sMap.codebook=sD.data(idx(1:munits));
