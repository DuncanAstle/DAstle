% GHSOM_DEMO1  Run a demo using a simple 2-dimensional dataset.
%
% Type 'ghsom_demo1' to start it.
%

% Copyright (c) 2001 by Elias Pampalk

% Version 1.0 Elias Pampalk 17072001
function ghsom_demo1

sigma = 0.1;
x=10;
y=10;
    
D = [reshape(repmat(1:x,y,1),1,x*y); repmat(y:-1:1,1,x)]'; % make 2-d rectangular grid
labels = strcat(int2str(D(:,1)),'|',int2str(D(:,2))); % label by grid position
D = D + randn(x*y,2)*sigma; % add a little noise

sD = som_data_struct(D,'name','Numbers 2D','labels',labels);

% comment: using extreme depth and breadth values to show orientation effects
ghM = ghsom_train(sD,'breadth',0.8,'depth',0.01,'tracking',0,'sub_layer_init','mirror'); 
ghM = ghsom_datalabels(ghM, sD);

f=figure; % plot dataset
set(f,'numbertitle','off');
set(f,'name','Demo 1: Testdata - Rectangular Grid with Noise');
plot(sD.data(:,1),sD.data(:,2),'.'); hold on
for i=1:x
    idx=((i-1)*y+1):((i-1)*y+x);
    plot(sD.data(idx,1),sD.data(idx,2),':'); 
end
for i=1:y
    idx=(1:x:x*y)+i-1;
    plot(sD.data(idx,1),sD.data(idx,2),':');
end
xlabel('X Value')
ylabel('Y Value')
title(['Testset Numbers 2D',' (\sigma=',num2str(sigma),')'])
axis([0 x+1 0 y+1])

% plot datalabels + grid and layer colors
f=figure; % plot dataset
set(f,'numbertitle','off');
set(f,'name','Demo 1: Data on Map (x|y - labels)');
ghV = ghsom_visualize_grid(ghM, 'layer');
ghsom_visualize_labels(ghM, ghV);

% plot component planes
f=figure; % plot dataset
set(f,'numbertitle','off');
set(f,'name','Demo 1: Component Planes');
ghsom_visualize_grid(ghM, 'component', [1:2]);
