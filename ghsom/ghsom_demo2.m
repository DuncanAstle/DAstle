% GHSOM_DEMO2  Run a demo using a simple 3-dimensional dataset.
%
% Type 'ghsom_demo2' to start it.
%
% The demo uses a three dimensional toy-set which can be used to show
% orientation effects as well as some hierarchy effects.
%

% Copyright (c) 2001 by Elias Pampalk

% Version 1.0 Elias Pampalk 17072001

function ghsom_demo2

% create dataset
sigma = 0.1;
x=10;
y=10;
D = [reshape(repmat(1:x,y,1),1,x*y); repmat(y:-1:1,1,x)]';
labels = strcat(int2str(D(:,1)),'|',int2str(D(:,2)));
D = D + randn(x*y,2)*sigma;
signes=sign((D(:,1)-mean(D(:,1))).^2-(D(:,2)-mean(D(:,2))).^2);
D(:,3)=signes.*sqrt(signes.*((D(:,1)-mean(D(:,1))).^2-(D(:,2)-mean(D(:,2))).^2));
sD=som_data_struct(D,'name','Numbers 3D','labels',labels);

% create GHSOM
% comment: using extreme depth and breadth values to show orientation effects
ghM = ghsom_train(sD,'breadth',0.8,'depth',0.1,'tracking',1,'sub_layer_init','dittenbach');
ghM = ghsom_datalabels(ghM, sD);

% visualize...

f=figure; % plot dataset
set(f,'numbertitle','off');
set(f,'name','Demo 2: Testdata - Rectangular Grid with Noise + Z Dimension');
    for i=1:y
        idx=((i-1)*y+1):((i-1)*y+x);
        plot3(sD.data(idx,1),sD.data(idx,2),sD.data(idx,3),'.'); hold on
        plot3(sD.data(idx,1),sD.data(idx,2),sD.data(idx,3),':'); 
    end
    for i=1:x
        idx=(1:x:x*y)+i-1;
        plot3(sD.data(idx,1),sD.data(idx,2),sD.data(idx,3),'.'); hold on
        plot3(sD.data(idx,1),sD.data(idx,2),sD.data(idx,3),':'); 
    end
    xlabel('X Value')
    ylabel('Y Value')
    zlabel('Z Value')
    title('Testset Numbers 3D')
    axis([0 x+1 0 y+1])

% plot datalabels + grid and layer colors
f=figure; % plot dataset
set(f,'numbertitle','off');
set(f,'name','Demo 2: Data on Map (x|y - labels)');
ghV = ghsom_visualize_grid(ghM, 'layer');
ghsom_visualize_labels(ghM, ghV);

% plot component planes
f=figure; % plot dataset
set(f,'numbertitle','off');
set(f,'name','Demo 2: Component Planes');
ghsom_visualize_grid(ghM, 'component', [1:7]);
