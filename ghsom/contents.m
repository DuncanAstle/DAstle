echo on
% GHSOM Toolbox for Matlab 5
% Extension of SOM Toolbox.
% 
% Create, train, visualize GHSOM.
%
% To get started try 'ghsom_demo'.
%
% FUNCTION OVERVIEW:
%
% Training:
% 
%   GHSOM_TRAIN                    all the training
%   GHSOM_TRAIN_GROW               only width growing
%   GHSOM_TRAIN_GROW_ORIENTATION   orientations of sublayers
%   GHSOM_TRAIN_GROW_SAMPLEINIT    sample initialization
%
% Labeling:
%
%   GHSOM_DATALABELS               label map using data labels
%   GHSOM_LABELSOM                 use LabelSOM method
%   GHSOM_WEBSOMLABELS             use WebSOM method
%
% Visualize:
%
%   GHSOM_VISUALIZE_GRID           creates 2D grid of GHSOM, can do component planes
%   GHSOM_VISUALIZE_LABELS         write labels into 2D grid
%   GHSOM_VISUALIZE_2DTS           visualize model vectors assuming they represent 2-dim time series
%
% External Interface:
%   
%   GHSOM_WRITE_CODS               SOM_PAK format output
%   GHSOM_READ_CODS                SOM_PAK format input
%
% Demonstrations:
%
%   GHSOM_DEMO                     use Kohonens animals
%   GHSOM_DEMO1                    use 2D numbers
%   GHSOM_DEMO2                    use 3D numbers
%   GHSOM_DEMO3                    use 3 normal distributed clusters (2D)s
%   GHSOM_DEMO4                    use music collection (77 pieces of music)
%   GHSOM_DEMO5                    use music performance data (2-dim time series)
%
% Other Stuff:
%
%   GHSOM_IDX2DATA                 get sData belonging to (sub)map knowing the indexes (ghMap.dataitems)
%   animals.dat                    Needed for GHSOM_DEMO
%   music.mat                      Needed for GHSOM_DEMO4
%   worms.mat                      Needed for GHSOM_DEMO5
%   GHSOM_MAP_DATA                 Use when loading ghMap from SOM_PAK format

% Copyright (c) 2001-2002 by Elias Pampalk, Alvin Chan

% Version 1.0 Alvin Chan 30042001
% Version 2.0 Elias Pampalk 17072001
% Version 2.1 Elias Pampalk 20052002
echo off
