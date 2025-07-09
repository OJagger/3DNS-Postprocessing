function [x,y,pitch] = periodic_boundary(blk)

[bx, by] = domain_boundary(blk);

if ispolycw(bx,by)
    bx = flip(bx);
    by = flip(by);
end

bx = reshape(bx,[],1);
by = reshape(by,[],1);

b = [bx by];
b = unique(b, 'rows', 'stable');

bx = b(:,1)';
by = b(:,2)';

xin = min(bx);
yin = by(bx == xin);
xout = max(bx);
yout = by(bx == xout);

iu = find(bx == xin & by == max(yin));

n = length(bx);
bx = [bx bx];
by = [by by];

bx = bx(iu:iu+n);
by = by(iu:iu+n);

il = find(bx == xin & by == min(yin));
ol = find(bx == xout & by == min(yout));

x = bx(il:ol);
y = by(il:ol);
pitch = by(1) - by(il);

end