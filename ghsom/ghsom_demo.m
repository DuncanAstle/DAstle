% GHSOM_DEMO  Run a demo using a simple toy dataset of animals.
%
% Type 'ghsom_demo' to start it. It returns the trained GHSOM.
%
% The demo uses a toy dataset published by Kohonen.
%
% See also GHSOM_DEMO1, GHSOM_DEMO2

% Copyright (c) 2001 by Elias Pampalk

% Version 1.0 Elias Pampalk 17072001

function ghMap = ghsom_demo

sD = som_read_data('animals.dat');

ghMap = ghsom_train(sD,'breadth',0.5,'depth',0.0001,'depth_dataitems',4,'layer_1st_init','lininit');
ghMap = ghsom_datalabels(ghMap, sD);

f=figure;
set(f,'numbertitle','off');
set(f,'name','Demo - Animals: GHSOM Mapping');
ghVisu = ghsom_visualize_grid(ghMap, 'layer');
ghsom_visualize_labels(ghMap, ghVisu);

f=figure;
set(f,'numbertitle','off');
set(f,'name','Demo - Animals: Component Planes');
ghsom_visualize_grid(ghMap, 'component', [1:13]);
