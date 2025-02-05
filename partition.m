function [pmin,bmin]=partition(dir,ncards)

bijk = load(fullfile(dir,'blockdims.txt'));
%ncards = 4

[nb,~] = size(bijk)
    
block_points = bijk(:,1).*bijk(:,2).*bijk(:,3)

total_points = sum(block_points)

points_per_card = total_points/ncards

block_array = perms([1:nb]);

%[B,I] = sort(block_points,'descend');

[nperms,~]=size(block_array);

err_min = 1e24;

for mm=1:5 % randomize
for i=1:nperms

num_blocks_left = nb;
mnow = 0;
p=[];
b=[];
nbb=[];
for n=1:ncards

% random number of blocks per card
ncards_left = ncards - n + 1;
nmax = num_blocks_left - ncards_left + 1;
if(n==ncards)
    num_blocks_now = num_blocks_left;
else
    num_blocks_now = randi(nmax);
end

num_blocks_left = num_blocks_left - num_blocks_now;


points = 0;    
for m=1:num_blocks_now
mnow = mnow + 1;    
nbnow = block_array(i,mnow);
points = points + block_points(nbnow);
p(n) = points;
b{n}(m)=nbnow;
nbb{n}=num_blocks_now;
end

end

err_now = sqrt(sum( (p - points_per_card).^2 ));
if(err_now<err_min)
pmin = p;
bmin = b;
nbbmin = nb;
err_min = err_now;
end

end
end

%[~,imin]=min(err);










