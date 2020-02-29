function ghsom_demo4
% GHSOM_DEMO4  Run a demo using a data set of 77 pieces of music.
%              See http://www.oefai.at/~elias/music for details.
%
% Type 'ghsom_demo4' to start it.
%
% This demo shows an application of the GHSOM to music collections.
% The obtained hierarchical organization is an intuitive interface to
% electronic music archives.
%

% Copyright (c) 2002 by Elias Pampalk

% Version 1.0 Elias Pampalk 20052002
%load worms
load music

for i=1:1
figure
ghMap =  ghsom_train(sData,'breadth',0.6,'depth',0.001,'tracking',0)
ghMap = ghsom_datalabels(ghMap, sData);
ghVisu = ghsom_visualize_grid(ghMap,'layer');
ghsom_visualize_labels(ghMap, ghVisu);
end