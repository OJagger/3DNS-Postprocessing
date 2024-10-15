function plot_mesh(blk, skip, f,fmt,label_blocks)
if nargin < 2 || isempty(skip)
    skip=1;
end
if nargin < 4 || isempty(fmt)
    fmt = 'k';
end
NB = length(blk.x);
if nargin < 3 || isempty(f)
    figure;
end
if nargin < 5 || isempty(label_blocks)
    label_blocks = false;
end
hold on
for i=1:NB
    
    xnew = blk.x{i};
    ynew = blk.y{i};
    [ni,nj] = size(xnew);
    is = round(linspace(1, ni, ni/skip));
    % is = 1:skip:ni;
    is = [is(1:end-1) ni];
    js = round(linspace(1, nj, nj/skip));
    % js = 1:skip:nj;
    % js = [js(1:end-1) nj];
    xnew = xnew(is, js);
    ynew = ynew(is,js);
    plot(xnew,ynew,fmt);
    plot(xnew',ynew',fmt);
    plot(blk.x{i}(1,:),blk.y{i}(1,:),'r');
    plot(blk.x{i}(end,:),blk.y{i}(end,:),'r');
    plot(blk.x{i}(:,1),blk.y{i}(:,1),'r');
    plot(blk.x{i}(:,end),blk.y{i}(:,end),'r');

end
if label_blocks
    for ib = 1:NB
        i = floor(blk.blockdims(ib,1)/2);
        j = floor(blk.blockdims(ib,2)/2);
        text(blk.x{ib}(i,j), blk.y{ib}(i,j), num2str(ib),"Color",'b');
    end
end
axis equal
end