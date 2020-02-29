% GHSOM_DEMO3  Run a demo comparing GHSOM and SOM using a simple 2-dimensional dataset.
%
% Type 'ghsom_demo3' to start it.
%
% The demo uses a three dimensional toy-set which can be used to show
% orientation effects as well as some hierarchy effects. It also tries to
% visualize magnification effects.
% There are 4 clusters. 2 bigger ones with a bigger distance to each other,
% and 2 smaller ones with a smaller distance and variance. It is interesting
% to see how these are represented with the GHSOM at different granularities.
%

% Copyright (c) 2001 by Elias Pampalk

% Version 1.0 Elias Pampalk 17072001

function ghsom_demo3

xy=zeros(2000,3);

xy(1:1000,1:2)=randn(1000,2)*.5+repmat([1 1],1000,1);
xy(1:1000,3)=eps;
xy(1001:2000,1:2)=randn(1000,2)*.2;
xy(1001:2000,3)=0;
xy(2001:2250,1:2)=randn(250,2)*.05+repmat([1.5, 1.5],250,1);
xy(2251:2500,1:2)=randn(250,2)*.05+repmat([1.8, 1.8],250,1);
xy(2001:2500,3)=0;

labels(1:1000)=repmat({'o'},1,1000);
labels(1001:2000)=repmat({'+'},1,1000);
labels(2001:2250)=repmat({'*'},1,250);
labels(2251:2500)=repmat({'x'},1,250);

sData = som_data_struct(xy,'labels',labels','comp_names',{'X','Y','Class'});

figure;
plot(xy(1:1000,1),xy(1:1000,2),'ok','markersize',3); hold on
plot(xy(1001:2000,1),xy(1001:2000,2),'+r','markersize',3);
plot(xy(2001:2250,1),xy(2001:2250,2),'*r','markersize',3);
plot(xy(2251:2500,1),xy(2251:2500,2),'xr','markersize',3);
xlabel('X');
ylabel('Y');

sMap=som_make(sData,'msize',[12 7],'tracking',0);
figure;
som_show(sMap,'comp','all','bar','none','edge','on','footnote','','subplots',[3,1]);
colormap gray

figure;
sMap = som_label(sMap,'clear','all');
sMap = som_autolabel(sMap,sData,'vote');
som_show(sMap,'footnote','','color',sMap.codebook(:,3),'bar','none','edge','on');
som_show_add('label',sMap);
colormap(repmat(linspace(0.5,1,64)',1,3))
shading faceted % get grid back

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ghMap = ghsom_train(sData,'breadth',0.5,'depth',0.007,'tracking',0,'sub_layer_init','dittenbach');
f=figure;
set(f,'numbertitle','off')
set(f,'name','Mirror')
ghVisu = ghsom_visualize_grid(ghMap,'component',[1:3]);
figure;
ghVisu = ghsom_visualize_grid(ghMap,'layer');
ghMap = ghsom_datalabels(ghMap, sData,'vote');
ghsom_visualize_labels(ghMap, ghVisu);