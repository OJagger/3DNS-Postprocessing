clear;

gridfldr = '/mnt/Tank/DNSData/Impeller/scratchStart2/300M/fine_grid/';%/mnt/Tank/DNSData/Impeller/scratchStart2/400M-2/2/grid2/';%runSetup/';%setUpArcher2/';%runSetupRfndMsh/att4/setUp/'; 
namePath = '../Data/turbin_rad-fine-TEST.dat';%Fine.dat';
gridRead = 'bin';

nmBlc = 2;
DIM = importdata([gridfldr,'blockdims.txt']);

grid = cell(2,1);

%I = 400;
iStart = 0;
iEnd = 300;%550; %420;
I = iEnd-iStart;
L = 0.009533835077005;

%% READ IN GRID AND BLOCKS

for nb = 1:nmBlc   
    
    ni = DIM(nb,1);
    nj = DIM(nb,2);
    nk = DIM(nb,3);
    
    
    if strcmp(gridRead,'txt')
        G = importdata([gridfldr,'grid_',num2str(nb),'.txt']);
        %G = reshape(G,[ni*nj*nk,3]); %reading griddata (vectorised)
        x = reshape(G(:,1),[ni,nj,nk]);
        y = reshape(G(:,2),[ni,nj,nk]);
        z = reshape(G(:,3),[ni,nj,nk]);
    elseif strcmp(gridRead,'bin')
        gridName = strcat(gridfldr,'grid_',num2str(nb));
        fid_g = fopen(gridName,'r'); %
        G = fread(fid_g,inf,'float64');       
        
        G = reshape(G,3,length(G)/3);  
        x = reshape(G(1,:),[ni,nj,nk]);
        y = reshape(G(2,:),[ni,nj,nk]);
        z = reshape(G(3,:),[ni,nj,nk]);
        fclose(fid_g); 
    else
        disp('Can\"t read grid')
        break
    end
    clear G
    
    grid{nb}.x = x;
    grid{nb}.y = y;
    grid{nb}.z = z;
    
    %[t,r,~] = cart2pol(z,y,x);

end

xb = cat(2,grid{2}.x,grid{1}.x(:,2:end,:));
yb = cat(2,grid{2}.y,grid{1}.y(:,2:end,:));
zb = cat(2,grid{2}.z,grid{1}.z(:,2:end,:));


dx = mean(mean(xb(2,:,:)-xb(1,:,:)))

xb = repmat(xb(1,:,:),I,1,1);
yb = repmat(yb(1,:,:),I,1,1);
zb = repmat(zb(1,:,:),I,1,1);

xfac = linspace(iStart,iEnd-1,I)';

%xfac = linspace(0,I-1,I)';

xfac = xfac*dx;

xfac = repmat(xfac,1,659,480); %479,360);

xb = xb+xfac;

clear grid x y z


figure(6)
plot3(xb(:,1,1),yb(:,1,1),zb(:,1,1))
hold on
plot3(xb(:,end,1),yb(:,end,1),zb(:,end,1))
plot3(xb(1,:,1),yb(1,:,1),zb(1,:,1))
plot3(xb(:,1,end),yb(:,1,end),zb(:,1,end))
hold on
plot3(xb(:,end,end),yb(:,end,end),zb(:,end,end))
plot3(xb(1,:,end),yb(1,:,end),zb(1,:,end))
plot3(squeeze(xb(1,1,:)),squeeze(yb(1,1,:)),squeeze(zb(1,1,:)))
hold on
plot3(squeeze(xb(end,end,:)),squeeze(yb(end,end,:)),squeeze(zb(end,end,:)))
plot3(squeeze(xb(1,end,:)),squeeze(yb(1,end,:)),squeeze(zb(1,end,:)))
plot3(squeeze(xb(end,1,:)),squeeze(yb(end,1,:)),squeeze(zb(end,1,:)))
plot3(squeeze(xb(end,:,1)),squeeze(yb(end,:,1)),squeeze(zb(end,:,1)))
plot3(squeeze(xb(end,:,end)),squeeze(yb(end,:,end)),squeeze(zb(end,:,end)))
xlabel('x')
ylabel('y')
zlabel('z')

x = linspace(min(xb(:)),L+min(xb(:)),100);
y = linspace(min(yb(:)),L+min(yb(:)),100);
z = linspace(min(zb(:)),L+min(zb(:)),100);
[X,Y,Z] = meshgrid(x,y,z);

plot3(X(:,1,1),Y(:,1,1),Z(:,1,1))
hold on
plot3(X(:,end,1),Y(:,end,1),Z(:,end,1))
plot3(X(1,:,1),Y(1,:,1),Z(1,:,1))
plot3(X(:,1,end),Y(:,1,end),Z(:,1,end))
plot3(X(:,end,end),Y(:,end,end),Z(:,end,end))
plot3(X(1,:,end),Y(1,:,end),Z(1,:,end))
plot3(squeeze(X(1,1,:)),squeeze(Y(1,1,:)),squeeze(Z(1,1,:)))
plot3(squeeze(X(end,end,:)),squeeze(Y(end,end,:)),squeeze(Z(end,end,:)))
plot3(squeeze(X(1,end,:)),squeeze(Y(1,end,:)),squeeze(Z(1,end,:)))
plot3(squeeze(X(end,1,:)),squeeze(Y(end,1,:)),squeeze(Z(end,1,:)))
plot3(squeeze(X(end,:,1)),squeeze(Y(end,:,1)),squeeze(Z(end,:,1)))
plot3(squeeze(X(end,:,end)),squeeze(Y(end,:,end)),squeeze(Z(end,:,end)))

clear X Y Z

[ni,nj,nk] = size(xb);

fid_f = fopen(namePath,'r'); %
A = fread(fid_f,inf,'float64');
A = reshape(A,3,nj,nk,ni);

u = squeeze(A(1,:,:,:));
v = squeeze(A(2,:,:,:));
w = squeeze(A(3,:,:,:));
u = permute(u,[3,1,2]);
v = permute(v,[3,1,2]);
w = permute(w,[3,1,2]);

[t,r] = cart2pol(zb,yb);

vr = cos(t).*w+sin(t).*v;
vt = cos(t).*v-sin(t).*w;

%pcolor(squeeze(u(1,:,:)))
%shading flat

figure(6)
surf(squeeze(xb(1,:,:)),squeeze(yb(1,:,:)),squeeze(zb(1,:,:)),squeeze(u(1,:,:)));
hold on
surf(squeeze(xb(end,:,:)),squeeze(yb(end,:,:)),squeeze(zb(end,:,:)),squeeze(u(end,:,:)));
surf(squeeze(xb(:,1,:)),squeeze(yb(:,1,:)),squeeze(zb(:,1,:)),squeeze(u(:,1,:)));
surf(squeeze(xb(:,end,:)),squeeze(yb(:,end,:)),squeeze(zb(:,end,:)),squeeze(u(:,end,:)));
surf(squeeze(xb(:,:,1)),squeeze(yb(:,:,1)),squeeze(zb(:,:,1)),squeeze(u(:,:,1)));
surf(squeeze(xb(:,:,end)),squeeze(yb(:,:,end)),squeeze(zb(:,:,end)),squeeze(u(:,:,end)));
colorbar()
%caxis([-25 25])
axis equal
shading flat
%{
figure(2)
surf(squeeze(xb(1,:,:)),squeeze(yb(1,:,:)),squeeze(zb(1,:,:)),squeeze(vt(1,:,:)));
hold on
surf(squeeze(xb(end,:,:)),squeeze(yb(end,:,:)),squeeze(zb(end,:,:)),squeeze(vt(end,:,:)));
surf(squeeze(xb(:,1,:)),squeeze(yb(:,1,:)),squeeze(zb(:,1,:)),squeeze(vt(:,1,:)));
surf(squeeze(xb(:,end,:)),squeeze(yb(:,end,:)),squeeze(zb(:,end,:)),squeeze(vt(:,end,:)));
surf(squeeze(xb(:,:,1)),squeeze(yb(:,:,1)),squeeze(zb(:,:,1)),squeeze(vt(:,:,1)));
surf(squeeze(xb(:,:,end)),squeeze(yb(:,:,end)),squeeze(zb(:,:,end)),squeeze(vt(:,:,end)));
colorbar()
caxis([-25 25])
axis equal
shading flat

figure(3)
surf(squeeze(xb(1,:,:)),squeeze(yb(1,:,:)),squeeze(zb(1,:,:)),squeeze(vr(1,:,:)));
hold on
surf(squeeze(xb(end,:,:)),squeeze(yb(end,:,:)),squeeze(zb(end,:,:)),squeeze(vr(end,:,:)));
surf(squeeze(xb(:,1,:)),squeeze(yb(:,1,:)),squeeze(zb(:,1,:)),squeeze(vr(:,1,:)));
surf(squeeze(xb(:,end,:)),squeeze(yb(:,end,:)),squeeze(zb(:,end,:)),squeeze(vr(:,end,:)));
surf(squeeze(xb(:,:,1)),squeeze(yb(:,:,1)),squeeze(zb(:,:,1)),squeeze(vr(:,:,1)));
surf(squeeze(xb(:,:,end)),squeeze(yb(:,:,end)),squeeze(zb(:,:,end)),squeeze(vr(:,:,end)));
colorbar()
caxis([-25 25])
axis equal
shading flat


figure(4)
surf(squeeze(xb(1,:,:)),squeeze(yb(1,:,:)),squeeze(zb(1,:,:)),squeeze(w(1,:,:)));
hold on
surf(squeeze(xb(end,:,:)),squeeze(yb(end,:,:)),squeeze(zb(end,:,:)),squeeze(w(end,:,:)));
surf(squeeze(xb(:,1,:)),squeeze(yb(:,1,:)),squeeze(zb(:,1,:)),squeeze(w(:,1,:)));
surf(squeeze(xb(:,end,:)),squeeze(yb(:,end,:)),squeeze(zb(:,end,:)),squeeze(w(:,end,:)));
surf(squeeze(xb(:,:,1)),squeeze(yb(:,:,1)),squeeze(zb(:,:,1)),squeeze(w(:,:,1)));
surf(squeeze(xb(:,:,end)),squeeze(yb(:,:,end)),squeeze(zb(:,:,end)),squeeze(w(:,:,end)));
colorbar()
caxis([-25 25])
axis equal
shading flat

figure()
scatter(u(:),vt(:))
axis equal
xlabel('$u$','interpreter','latex','fontsize',15)
ylabel('$v_{\theta}$','interpreter','latex','fontsize',15)

figure()
scatter(u(:),vr(:))
axis equal
xlabel('$u$','interpreter','latex','fontsize',15)
ylabel('$v_{r}$','interpreter','latex','fontsize',15)

figure()
scatter(vr(:),vt(:))
axis equal
xlabel('$v_{r}$','interpreter','latex','fontsize',15)
ylabel('$v_{\theta}$','interpreter','latex','fontsize',15)

figure()
scatter(u(:),v(:))
axis equal
xlabel('$u$','interpreter','latex','fontsize',15)
ylabel('$v$','interpreter','latex','fontsize',15)

figure()
scatter(u(:),w(:))
axis equal
xlabel('$u$','interpreter','latex','fontsize',15)
ylabel('$w$','interpreter','latex','fontsize',15)

figure()
scatter(v(:),w(:))
axis equal
xlabel('$v$','interpreter','latex','fontsize',15)
ylabel('$w$','interpreter','latex','fontsize',15)
axis equal
%}

%{
figure
ud = u(:,1,:)-u(:,end,:);
 imshow(squeeze(ud))
shading flat
colorbar
vrd =vr(:,1,:)-vr(:,end,:);
figure
imshow(squeeze(vrd))
vtd =vt(:,1,:)-vt(:,end,:);
figure
imshow(squeeze(vtd))

ud =u(1,:,:)-u(end,:,:);
figure
imshow(squeeze(ud))
%}
%{
[zb(1,1,1),zb(1,end,1),zb(1,1,end),zb(1,end,end)]
Z=zb(1,:,:);
Y=yb(1,:,:);
Z=squeeze(Z)
Y=squeeze(Y);
X=squeeze(xb(1,:,:));
yMax=max(Y(:))
yMin=min(Y(:))
zMax=max(Z(:))
zMin=min(Z(:))
xyz=cat(3,X,Y,Z);
dI=xyz(2:end,:,:)-xyz(1:end-1,:,:);
lI=vecnorm(dI,2,3);
lMinI=min(lI(:))
which(lI==0)
find(lI==0)
ind2sub([899,240],450)
xyz(450,1,:)-xyz(451,1,:)
dJ=xyz(:,2:end,:)-xyz(:,1:end-1,:);
lJ=vecnorm(dJ,2,3);
lMinJ=min(lJ(:))
xyz2=xyz([1:450,452:end],:,:);
dI=xyz2(2:end,:,:)-xyz2(1:end-1,:,:);
lI=vecnorm(dI,2,3);
lMinI=min(lI(:))
lMin=min(lMinI,lMinJ)
(yMax-yMin)*(zMax-zMin)/lMin^2
yl = linspace(yMin,yMax,900);
zl = linspace(zMin,zMax,900);
yMax-yMin
zMax-zMin
[t,r] = cart2pol(zb,yb);
min(min(min(r)))*2*pi/9
minL = min(min(min(r)))*2*pi/9
maxL = zMax-zMin
maxL/lMin
%}
%}