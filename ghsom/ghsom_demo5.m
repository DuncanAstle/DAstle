function ghsom_demo5
% GHSOM_DEMO5  Run a demo using a data set of multivariate time series sequences.
%              See http://www.oefai.at/music for details. The data represents
%              musical performances of classical piano pieces, in particular, of
%              Chopin. The 40 dimensions represent: 1-20: tempo [bpm], 21-40: 
%              loudness [sone].
%              See also: "Hierarchical Clustering and Structuring of Data with
%                         Self-Organizing Maps", E. Pampalk, G. Widmer, and A. Chan
%              for a detailed discussion of the results.
%
% Type 'ghsom_demo5' to start it.
%

% Copyright (c) 2002 by Elias Pampalk

% Version 1.0 Elias Pampalk 02062002

load worms

ghMap =  ghsom_train(sData,'breadth',0.7,'depth',0.05,'tracking',2,'depth_dataitems',10)

figure
ghVisu = ghsom_visualize_grid(ghMap,'layer',1);
ghsom_visualize_2dts(ghMap,ghVisu,sData,1);

figure
ghVisu = ghsom_visualize_grid(ghMap,'layer');
ghsom_visualize_2dts(ghMap,ghVisu,sData);
