% GHSOM_TRAIN_GROW_ORIENTATION  Initialize a (sub) map for the GHSOM with a certain orientation 
%
%  sMap = ghsom_train_grow_orientation(parent_codebooks, no_edge, method)
%
%  Input and output arguments: 
%   parent_codebooks (matrix) codebooks (model vectors) parent layer, size 9 x dim
%   no_edge          (struct) indicates if some neighbors are missing
%   method           (string) orientation type choose 'mirror' (default) or 'dittenbach'
%   sMap             (struct) SOM Toolbox map structure
%
%   
% Try 'type ghsom_orientation' for more details.
% See also GHSOM_TRAIN, GHSOM_TRAIN_GROW

%%%%%%%%%%%%% DETAILED DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ghsom_train_grow_orientation
% 
% PURPOSE
%
% Keep the orientation of the upper level maps on the lower level maps.
% Creates an initialization of a 2x2 SOM based on the parent unit and its neighbors on
% the parent map.
%
% DESCRIPTION
%
% Currently two methods are implemented: 'mirror' and 'dittenbach'. The second refers
% to the method published by Michael Dittenbach 2001 with the modification that
% the the formulas of the type: 'a4 = a + ((b-a)+(d-a)+(e-a))/3' have been changed to
% 'a4 = a + ((b-a)+(d-a)+(e-a))/6' calculated as 'a4 = mean([a,mean(b,d,e)])'.
% Dittenbachs published method uses 8 neighbors of the parent unit and if neighbors
% are missing (if the parent unit is located at some edge of the map) then the newly
% created 2x2 SOM will not have any map units outside of the area covered py the
% parent map. This is done for example with the formula of the type: 'a1=a', which
% initializes the 'corner' of the sub map to the value of the parent unit.
% The method named 'mirror' has two main differences: first of all only 4 neighbors
% of the parent unit are used. The idea is that the corner neigbhors are sqrt(2) times
% further away than the immidiate 4 neighbors and - in case of a finetuned (small) SOM
% might be quite far away in the input space. The values (in the notation of Dittenbach)
% are calculated as a4 = mean(a,b,d). The other main difference is how missing units
% are treated i.e. when the parent unit is alligned to some edge. Assuming that the
% SOM localy approximates a linear hyperplane the missing units are estimated by 
% 'mirroring' the existing ones. For example if the left unit is missing, then
% a virtual right unit is created reflecting the left one on the otherside of the
% center unit.
%
% REFERENCES
%
% M. Dittenbach, A. Rauber, and D. Merkl. , "Recent Advances with the Growing 
%    Hierarchical Self-Organizing Map", N. Allison, H. Yin, L. Allison, 
%    and J. Slack, editors, Advances in Self-Organising Maps, pp. 140-145, 
%    June 13-15, 2001, Lincoln, UK, Springer-Verlag.
%
% Copyright (c) 2001 by Elias Pampalk

% Version 1.1 Elias Pampalk 29102002 fixed some bug
% Version 1.0 Elias Pampalk 17072001

function sMap = ghsom_train_grow_orientation(parent_codebooks, no_edge, method)
% 1 4 7
% 2 5 8
% 3 6 9

switch method
case 'dittenbach',
    if no_edge.down & no_edge.up & no_edge.left & no_edge.right, % 1:9
        child_codebook(1,:)=mean([parent_codebooks(5,:);mean([parent_codebooks(1,:);parent_codebooks(2,:);parent_codebooks(4,:)])]); % upper left
        child_codebook(2,:)=mean([parent_codebooks(5,:);mean([parent_codebooks(2,:);parent_codebooks(3,:);parent_codebooks(6,:)])]); % lower left
        child_codebook(3,:)=mean([parent_codebooks(5,:);mean([parent_codebooks(4,:);parent_codebooks(7,:);parent_codebooks(8,:)])]); % upper right
        child_codebook(4,:)=mean([parent_codebooks(5,:);mean([parent_codebooks(6,:);parent_codebooks(8,:);parent_codebooks(9,:)])]); % lower right
    elseif ~no_edge.left & ~no_edge.up, % 5,6,8,9
        child_codebook(1,:)=parent_codebooks(5,:); % upper left
        child_codebook(2,:)=mean([parent_codebooks(5,:);parent_codebooks(6,:)]); % lower left
        child_codebook(3,:)=mean([parent_codebooks(5,:);parent_codebooks(8,:)]); % upper right
        child_codebook(4,:)=mean([parent_codebooks(5,:);mean([parent_codebooks(6,:);parent_codebooks(8,:);parent_codebooks(9,:)])]); % lower right
    elseif ~no_edge.left & no_edge.down & no_edge.up, %4:9
        child_codebook(1,:)=mean([parent_codebooks(5,:);parent_codebooks(4,:)]); % upper left
        child_codebook(2,:)=mean([parent_codebooks(5,:);parent_codebooks(6,:)]); % lower left
        child_codebook(3,:)=mean([parent_codebooks(5,:);mean([parent_codebooks(4,:);parent_codebooks(7,:);parent_codebooks(8,:)])]); % upper right
        child_codebook(4,:)=mean([parent_codebooks(5,:);mean([parent_codebooks(6,:);parent_codebooks(8,:);parent_codebooks(9,:)])]); % lower right        
    elseif ~no_edge.left & ~no_edge.down, % 4,5,7,8
        child_codebook(1,:)=mean([parent_codebooks(5,:);parent_codebooks(4,:)]); % upper left
        child_codebook(2,:)=parent_codebooks(5,:); % lower left
        child_codebook(4,:)=mean([parent_codebooks(5,:);parent_codebooks(8,:)]); % upper right
        child_codebook(3,:)=mean([parent_codebooks(5,:);mean([parent_codebooks(4,:);parent_codebooks(7,:);parent_codebooks(8,:)])]); % lower right
    elseif ~no_edge.up & no_edge.left & no_edge.right, % 2,3,5,6,8,9
        child_codebook(1,:)=mean([parent_codebooks(5,:);parent_codebooks(2,:)]); % upper left
        child_codebook(2,:)=mean([parent_codebooks(5,:);mean([parent_codebooks(2,:);parent_codebooks(3,:);parent_codebooks(6,:)])]);
        child_codebook(3,:)=mean([parent_codebooks(5,:);parent_codebooks(8,:)]); 
        child_codebook(4,:)=mean([parent_codebooks(5,:);mean([parent_codebooks(6,:);parent_codebooks(8,:);parent_codebooks(9,:)])]); % lower right        
    elseif ~no_edge.down & no_edge.left & no_edge.right, % 1,2,4,5,7,8
        child_codebook(1,:)=mean([parent_codebooks(5,:);mean([parent_codebooks(1,:);parent_codebooks(2,:);parent_codebooks(4,:)])]); % upper left
        child_codebook(2,:)=mean([parent_codebooks(5,:);parent_codebooks(2,:)]); % lower left
        child_codebook(3,:)=mean([parent_codebooks(5,:);mean([parent_codebooks(4,:);parent_codebooks(7,:);parent_codebooks(8,:)])]); % upper right
        child_codebook(4,:)=mean([parent_codebooks(5,:);parent_codebooks(8,:)]); % lower right        
    elseif ~no_edge.up & ~no_edge.right, % 2,3,5,6
        child_codebook(1,:)=mean([parent_codebooks(5,:);parent_codebooks(2,:)]); % upper left
        child_codebook(2,:)=mean([parent_codebooks(5,:);mean([parent_codebooks(2,:);parent_codebooks(3,:);parent_codebooks(6,:)])]); % upper right
        child_codebook(3,:)=parent_codebooks(5,:); % lower left
        child_codebook(4,:)=mean([parent_codebooks(5,:);parent_codebooks(6,:)]); % lower right        
    elseif ~no_edge.right & no_edge.down & no_edge.up, % 1:6
        child_codebook(1,:)=mean([parent_codebooks(5,:);mean([parent_codebooks(1,:);parent_codebooks(2,:);parent_codebooks(4,:)])]); % upper left
        child_codebook(2,:)=mean([parent_codebooks(5,:);mean([parent_codebooks(2,:);parent_codebooks(3,:);parent_codebooks(6,:)])]); % lower left
        child_codebook(3,:)=mean([parent_codebooks(5,:);parent_codebooks(4,:)]); % upper right
        child_codebook(4,:)=mean([parent_codebooks(5,:);parent_codebooks(6,:)]); % lower right        
    elseif ~no_edge.right & ~no_edge.down, % 1,2,4,5
        child_codebook(1,:)=mean([parent_codebooks(5,:);mean([parent_codebooks(1,:);parent_codebooks(2,:);parent_codebooks(4,:)])]); % upper left
        child_codebook(2,:)=mean([parent_codebooks(5,:);parent_codebooks(2,:)]); % lower left
        child_codebook(3,:)=mean([parent_codebooks(5,:);parent_codebooks(4,:)]); % upper right
        child_codebook(4,:)=parent_codebooks(5,:); % lower right
    end
otherwise, % 'mirror' which is default
    
    % replace missing sides. (max 2 missing, never 2 opposite ones are missing at the same time)    
    if ~no_edge.left % left
        parent_codebooks(2,:)=parent_codebooks(5,:)+parent_codebooks(5,:)-parent_codebooks(8,:);   
    end
    if ~no_edge.up % up-center
        parent_codebooks(4,:)=parent_codebooks(5,:)+parent_codebooks(5,:)-parent_codebooks(6,:);
    end
    if ~no_edge.down % down-center
        parent_codebooks(6,:)=parent_codebooks(5,:)+parent_codebooks(5,:)-parent_codebooks(4,:);
    end
    if ~no_edge.right % right
        parent_codebooks(8,:)=parent_codebooks(5,:)+parent_codebooks(5,:)-parent_codebooks(2,:);
    end
    
    % calculate initial codebooks as average of center + 2 closest neighbors
    child_codebook(1,:)=mean([parent_codebooks(5,:);parent_codebooks(2,:);parent_codebooks(4,:)]); % upper left
    child_codebook(2,:)=mean([parent_codebooks(5,:);parent_codebooks(2,:);parent_codebooks(6,:)]); % lower left
    child_codebook(3,:)=mean([parent_codebooks(5,:);parent_codebooks(4,:);parent_codebooks(8,:)]); % upper right
    child_codebook(4,:)=mean([parent_codebooks(5,:);parent_codebooks(6,:);parent_codebooks(8,:)]); % lower right
    
end

sMap = som_map_struct(size(parent_codebooks,2),'msize',[2 2],'sheet','rect');
sMap.codebook = child_codebook;
