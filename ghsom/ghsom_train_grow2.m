function sMap = ghsom_train_grow2(sData, parent_codebooks, no_edge)

NEIGH = 'gaussian';
%r = linspace(0.5,0.4,15);
r = linspace(0.7,0.2,15);

if nargin>1,
    sMap = ghsom_train_grow_orientation(parent_codebooks, no_edge,'mirror');
else
    %sMap = som_lininit(sData,'msize',[2 2],'rect','sheet');
    sMap = som_randinit(sData,'msize',[2 2],'rect','sheet');
end

% first training
sMap = som_batchtrain(sMap,sData,'tracking',0,'radius',r,'neigh',NEIGH);

% rot 90 codebook
%idx = rot90(reshape(1:4,2,2));
%sMap.codebook = sMap.codebook(idx(:),:);

error_old = REALMAX;
sMap_new = sMap;
k = 1;
while 1,
    [sMap_new_new, error_new] = grow(sMap_new,sData);
    disp([num2str(k),': ',num2str(error_new)]); k=k+1;
    if error_new < error_old,
        sMap = sMap_new;

        r = linspace(0.7+k-1,0.2,35+k*5);
        %r = linspace(0.5,0.4,10);
        sMap_new = som_batchtrain(sMap_new_new,sData,'tracking',0,'radius',r);
        error_old = error_new;
    else
        break;
    end
end

function [sMap_new, error] = grow(sMap,sData)
% calculate the number of items which would be mapped to intermediate virtual units, inkl border units

munits = size(sMap.codebook,1);
dims = size(sMap.codebook,2);
rows = sMap.topol.msize(1);
cols = sMap.topol.msize(2);

[bmus qerr] = som_bmus(sMap, sData); % which data subset belongs to each unit
unit_qerr = zeros(munits,1);
unit_mapped = zeros(munits,1);

for i = 1:munits,
    unit_qerr(i) = sum(qerr(find(bmus==i)));
    unit_mapped(i) = sum(bmus==i);
end
idx = find(unit_mapped>0);
unit_qerr(idx) = unit_qerr(idx) ./ unit_mapped(idx);

neighs = som_unit_neighs(sMap.topol);

M = sMap.codebook;

m5 = zeros(munits,5,dims); % model vecs + virtual vecs (center, up, down, left, right)
b5 = zeros(munits,5); % bmus count for m5 (center, up, down, left, right)

m5(:,1,:) = M;

for i=1:munits,
    if sum(neighs(i,:)) == 2, % corner
        if i == 1, % need up and left
            m5(i,2,:) = .5*(3*M(i,:)-M(i+1,:)); % virtual up
            m5(i,3,:) = .5*(M(i,:)+M(i+1,:)); % down
            m5(i,4,:) = .5*(3*M(i,:)-M(i+rows,:)); % virtual left
            m5(i,5,:) = .5*(M(i,:)+M(i+rows,:)); % right
        elseif i == rows, % need down and left
            m5(i,2,:) = .5*(M(i,:)+M(i-1,:)); % up
            m5(i,3,:) = .5*(3*M(i,:)-M(i-1,:)); % virtual down
            m5(i,4,:) = .5*(3*M(i,:)-M(i+rows,:)); % virtual left
            m5(i,5,:) = .5*(M(i,:)+M(i+rows,:)); % right
        elseif i == (cols-1)*rows+1, % need up and right
            m5(i,2,:) = .5*(3*M(i,:)-M(i+1,:)); % virtual up
            m5(i,3,:) = .5*(M(i,:)+M(i+1,:)); % down
            m5(i,4,:) = .5*(M(i,:)+M(i-rows,:)); % left
            m5(i,5,:) = .5*(3*M(i,:)-M(i-rows,:)); % virtual right
        else, % need down and right
            m5(i,2,:) = .5*(M(i,:)+M(i-1,:)); % up
            m5(i,3,:) = .5*(3*M(i,:)-M(i-1,:)); % virtual down
            m5(i,4,:) = .5*(M(i,:)+M(i-rows,:)); % left
            m5(i,5,:) = .5*(3*M(i,:)-M(i-rows,:)); % virtual right
        end
    elseif sum(neighs(i,:)) == 3, % border
        if mod(i,rows)==1, % need up
            m5(i,2,:) = .5*(3*M(i,:)-M(i+1,:)); % virtual up
            m5(i,3,:) = .5*(M(i,:)+M(i+1,:)); % down
            m5(i,4,:) = .5*(M(i,:)+M(i-rows,:)); % left
            m5(i,5,:) = .5*(M(i,:)+M(i+rows,:)); % right
        elseif mod(i,rows)==0, % need down
            m5(i,2,:) = .5*(M(i,:)+M(i-1,:)); % up
            m5(i,3,:) = .5*(3*M(i,:)-M(i-1,:)); % virtual down
            m5(i,4,:) = .5*(M(i,:)+M(i-rows,:)); % left
            m5(i,5,:) = .5*(M(i,:)+M(i+rows,:)); % right
        elseif i<rows, % need left
            m5(i,2,:) = .5*(M(i,:)+M(i-1,:)); % up
            m5(i,3,:) = .5*(M(i,:)+M(i+1,:)); % down
            m5(i,4,:) = .5*(3*M(i,:)-M(i+rows,:)); % virtual left
            m5(i,5,:) = .5*(M(i,:)+M(i+rows,:)); % right
        else, % need right
            m5(i,2,:) = .5*(M(i,:)+M(i-1,:)); % up
            m5(i,3,:) = .5*(M(i,:)+M(i+1,:)); % down
            m5(i,4,:) = .5*(M(i,:)+M(i-rows,:)); % left
            m5(i,5,:) = .5*(3*M(i,:)-M(i-rows,:)); % virtual right
        end
    else, % normal
        m5(i,2,:) = .5*(M(i,:)+M(i-1,:)); % up
        m5(i,3,:) = .5*(M(i,:)+M(i+1,:)); % down
        m5(i,4,:) = .5*(M(i,:)+M(i-rows,:)); % left
        m5(i,5,:) = .5*(M(i,:)+M(i+rows,:)); % right
    end
end

for i=1:munits,
    subdata = sData.data(bmus==i,:);
    bmuss = som_bmus(squeeze(m5(i,:,:)),subdata);
    for j=1:5,
        b5(i,j) = sum(bmuss == j);
    end
end

error = sum(sum(b5(:,2:5)))/sum(b5(:,1));

% limit TENSION
tr = zeros(rows-1,1);
tc = zeros(cols-1,1);
for i=0:rows-2,
    tr(i+1) = mean(sqrt(sum((M((1:rows:end)+i,:)-M((1:rows:end)+i+1,:)).^2,2)));
end
for i=0:cols-2,
    tc(i+1) = mean(sqrt(sum(M((1:rows)+rows*i,:)-M((1:rows)+rows*(i+1),:),2).^2));
end

%mtr = mean(tr)/mean(tc);
%mtc = mean(tc)/mean(tr);
%mtr = median(tr)/median(tc);
%mtc = median(tc)/median(tr);
mtr = tr/median(tc);
mtc = tc/median(tr);

mtr = [mtr(1); mtr; mtr(end)];
mtc = [mtc(1); mtc; mtc(end)];

mxtr = max(tr);
tr = tr/max(tc);
tc = tc/mxtr;

% insert_idx (best row or col)
r = zeros(rows+1,1);
r(1) = sum(b5(1:rows:end,2));
r(rows+1) = sum(b5(rows:rows:end,3));
c = zeros(cols+1,1);
c(1) = sum(b5(1:rows,4));
c(cols+1) = sum(b5((cols-1)*rows+1:end,5));

for i=0:rows-2,
    r(i+2) = sum(b5((1:rows:end)+i,3))+sum(b5((1:rows:end)+i+1,2));
end

for i=0:cols-2,
    c(i+2) = sum(b5((1:rows)+rows*i,5))+sum(b5((1:rows)+rows*(i+1),4));
end

r = r./cols.*mtr;
c = c./rows.*mtc;

if max(r)>max(c), % mean(r)>mean(c) !?
    idx = find(r==max(r));
    insert_idx=[idx(1),0];
else
    idx = find(c==max(c));
    insert_idx=[0,idx(1)];
end

sMap_new = sMap;
if insert_idx(1), % insert row
    sMap_new.codebook = zeros((rows+1)*cols,dims);
    idx = ones(rows+1,cols);
    sMap_new.topol.msize(1)=rows+1;
    if insert_idx(1)==1,
        sMap_new.codebook(1:rows+1:end,:) = squeeze(m5(1:rows:end,2,:));
        idx(1:rows+1:end) = 0;
    else,
        sMap_new.codebook((1:rows+1:end)+insert_idx(1)-1,:) = squeeze(m5((1:rows:end)+(insert_idx(1)-2),3,:));
        idx((1:rows+1:end)+insert_idx(1)-1) = 0;
    end
else % insert column
    sMap_new.codebook = zeros((cols+1)*rows,dims);
    idx = ones(rows,cols+1);
    sMap_new.topol.msize(2)=cols+1;
    if insert_idx(2)==1,
        sMap_new.codebook(1:rows,:) = squeeze(m5(1:rows,4,:));
        idx(1:rows) = 0;
    else,
        sMap_new.codebook((1:rows)+(insert_idx(2)-1)*rows,:) = squeeze(m5((1:rows)+(insert_idx(2)-2)*rows,5,:));
        idx((1:rows)+(insert_idx(2)-1)*rows) = 0;
    end
end
sMap_new.codebook(logical(idx(:)),:) = sMap.codebook;
sMap_new.labels = cell(size(sMap_new.codebook,1),1); % Add extra label fields for new units