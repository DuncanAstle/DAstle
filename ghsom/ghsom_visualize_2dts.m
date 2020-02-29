% GHSOM_VISUALIZE_2DTS  Plot 2-Dimensional Time Series on top of grid. 
%
%  ghVisu = ghsom_visualize_2dts(ghMap, ghVisu, sData, layer)
%
%  ghMap = ghsom_train(Data);
%  ghVisu = ghsom_visualize_grid(ghMap,'layer',2);
%  ghsom_visualize_2dts(ghMap,ghVisu,sData,2);

% Copyright (c) 2002 by Elias Pampalk, Alvin Chan

% Version 1.0 Elias Pampalk 30052002

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ghsom_visualize_2dts(ghMap, ghVisu, sData, layer)

if nargin<4,
    layer=0;
end

mini = min(min(sData.data)); maxi = max(max(sData.data));

min_scale_x = REALMAX;
min_scale_y = REALMAX;
matched = [];
k=1;
for map_idx=1:length(ghMap.sMap),
    if ~layer | ghMap.layer(map_idx)<=layer,
        bmus = som_bmus(ghMap.sMap{map_idx},sData.data(ghMap.dataitems{map_idx},:));
        if ghVisu.unit_size_xy{map_idx}(1) < min_scale_x,
            min_scale_x=ghVisu.unit_size_xy{map_idx}(1);
        end
        if ghVisu.unit_size_xy{map_idx}(2) < min_scale_y,
            min_scale_y=ghVisu.unit_size_xy{map_idx}(2);
        end
        for unit_idx=1:prod(ghMap.sMap{map_idx}.topol.msize),
            if (~layer & no_sub_unit(map_idx,unit_idx,ghMap.parent,ghMap.parent_unit)) | (layer & ghMap.layer(map_idx)==layer), 
                matched(k) = size(sData.data(ghMap.dataitems{map_idx}(bmus==unit_idx),:),1);
                k=k+1;
            end
        end
    end
end

for map_idx=1:length(ghMap.sMap),
    bmus = som_bmus(ghMap.sMap{map_idx},sData.data(ghMap.dataitems{map_idx},:));
    for unit_idx=1:prod(ghMap.sMap{map_idx}.topol.msize),
        if (~layer & no_sub_unit(map_idx,unit_idx,ghMap.parent,ghMap.parent_unit)) | (layer & ghMap.layer(map_idx)==layer), 
            visu_mapunit(...
                ghVisu.coordinates_xy{map_idx}(unit_idx,1),...
                ghVisu.coordinates_xy{map_idx}(unit_idx,2),...
                min_scale_x, ...
                min_scale_y, ...
                ... ghVisu.unit_size_xy{map_idx}(1),...
                ... ghVisu.unit_size_xy{map_idx}(2),...
                ghMap.sMap{map_idx}.codebook(unit_idx,:),...
                sData.data(ghMap.dataitems{map_idx}(bmus==unit_idx),:), mini,maxi,mean(matched));
        end
    end    
end

function is_true = no_sub_unit(map_idx, unit_idx, parentmap, parent_unit)
is_true=1;
children=find(parentmap==map_idx);
if ~isempty(children),
    if ~isempty(find(parent_unit(children)==unit_idx)),
        is_true=0;
    end
end

function visu_mapunit(pos_x, pos_y, scale_x, scale_y, c, D, mini, maxi,mean_matched)
matched = size(D,1);
if matched > 0,
    
    b = (mini+maxi)/2;
    s = maxi-mini;
    
    %scale = min(scale_x, scale_y)*0.85;
    scale_x=scale_x*0.85;
    scale_y=scale_y*0.85;
    
    px_mini = (mini-b)/s*scale_x + pos_x;
    py_mini = (mini-b)/s*scale_y + pos_y;
    px_maxi = (maxi-b)/s*scale_x + pos_x;
    py_maxi = (maxi-b)/s*scale_y + pos_y;
    
    BIG_PATCH_RES=200;
    big_xx = linspace(0,2*pi,BIG_PATCH_RES); % patch resolution
    big_PX = sin(big_xx); big_PY = cos(big_xx);
    
    h=patch(big_PX*scale_x/2*1.1 + pos_x,big_PY*scale_y/2*1.1 + pos_y,[1 1 1]);
    set(h,'linestyle','none')
    
    %patch([px_mini px_mini; px_mini px_maxi; px_maxi px_maxi; px_maxi px_mini], [py_mini py_maxi; py_maxi py_maxi; py_maxi py_mini; py_mini py_mini],[1 1 1]);
    set(h,'linestyle','none')
    
    zx = (0-b)/s*scale_x + pos_x;
    zy = (0-b)/s*scale_y + pos_y;
    
    %plot([px_mini px_maxi],[zy zy],':','color',[.9 .9 .9]); hold on
    %plot([zx zx], [py_mini py_maxi],':','color',[.9 .9 .9])
    
    px_c = (c(1:end/2)-b)/s*scale_x + pos_x;
    py_c = (c(end/2+1:end)-b)/s*scale_y + pos_y;
    
    Q=.5; % quantile
    sseg = D-repmat(c,matched,1);
    dist = sqrt(sseg(:,1:end/2).^2+sseg(:,end/2+1:end).^2);
    [dist dist_idx] = sort(dist,1);
    q = dist(max(round(matched*Q),1),:);
    
    PATCH_RES=20;
    xx = linspace(0,2*pi,PATCH_RES); % patch resolution
    PX = sin(xx); PY = cos(xx);
    
    PX2 = repmat(PX,length(c)/2,1)'.*repmat(q,PATCH_RES,1) + repmat(c(1:end/2),PATCH_RES,1);
    PY2 = repmat(PY,length(c)/2,1)'.*repmat(q,PATCH_RES,1) + repmat(c(end/2+1:end),PATCH_RES,1);
    
    PX2 = (PX2-b)/s*scale_x + pos_x;
    PY2 = (PY2-b)/s*scale_y + pos_y;
    
    p = patch(PX2,PY2,[0.9 0.9 1]); hold on
    set(p,'linestyle','none')

    plot(px_c,py_c,'-b'); hold on % worm colors
    plot(px_c(end),py_c(end),'.r')

    h=line(big_PX*scale_x/2*1.1 + pos_x,big_PY*scale_y/2*1.1 + pos_y);
    %h=line([px_mini px_mini; px_mini px_maxi; px_maxi px_maxi; px_maxi px_mini], [py_mini py_maxi; py_maxi py_maxi; py_maxi py_mini; py_mini py_mini]);
    set(h,'color',[.7 .7 .7])
    
    if 0,
        if matched>mean_matched*.66,
            h=line(big_PX*scale_x/2*1.1 + pos_x,big_PY*scale_y/2*1.1 + pos_y);
            set(h,'color',[.7 .7 .7])
        end
        if matched>mean_matched*1.33,
            h=line(big_PX*scale_x/2*1.2 + pos_x,big_PY*scale_y/2*1.2 + pos_y);
            set(h,'color',[.7 .7 .7])
        end  
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%