function mesh_analysis(blk, skip, Re)

if nargin < 3, Re=800e3; end

if nargin<2, skip=8; end

temp_slash = '/'; if ispc, temp_slash = '\'; end

C = colororder;

NB = length(blk.x);
try
    npp = blk.npp;
catch
    npp = 1;
end

nij_pts = 0;
figure
hold on
for i=1:NB
    x = blk.x{i};
    y = blk.y{i};
    x = x(:,:,1);
    y = y(:,:,1);
    procdims(i,:) = size(x)/npp;
    ni = size(x,1);
    nj = size(x,2);
    nij_pts = nij_pts + ni*nj;
    nk = blk.nk;
    for j=[1:skip:ni ni]
        if (j==1) || (j==ni)
            plot(x(j,:),y(j,:),'r','LineWidth',1)
        else
            plot(x(j,:),y(j,:),'k')
        end
    end
    for j=[1:skip:nj nj]
        if (j==1) || (j==nj)
            plot(x(:,j),y(:,j),'r','LineWidth',1)
        else
            plot(x(:,j),y(:,j),'k')
        end
    end
%     plot(x(1,:),y(1,:),'r','LineWidth',1)
%     plot(x(end,:),y(end,:),'r','LineWidth',1)
%     plot(x(:,1),y(:,1),'r','LineWidth',1)
%     plot(x(:,end),y(:,end),'r','LineWidth',1)
end

axis equal
%%

% data = readmatrix([case_name temp_slash 'probe.txt']);
% for i=1:size(data,1)
%     probe{i}.nb = data(i,1);
%     probe{i}.ni = data(i,2);
%     probe{i}.nj = data(i,3);
%     probe{i}.nk = data(i,4);
% end
% 
% 
% for i=1:length(probe)
%     probe_plot(i) = scatter(blk{probe{i}.nb}.x(probe{i}.ni, probe{i}.nj),blk{probe{i}.nb}.y(probe{i}.ni, probe{i}.nj));
%     probe_legend(i) = string(sprintf('Probe %d',i));
% end
% 
% trip.nb = 6;
% trip.i = 20;
% trip.j = size(blk{trip.nb}.x,2);
% trip.x = blk{probe{i}.nb}.x(probe{i}.ni, probe{i}.nj);
% trip.y = blk{probe{i}.nb}.y(probe{i}.ni, probe{i}.nj);
% trip
% 
% probe_plot(length(probe)+1) = scatter(blk{trip.nb}.x(trip.i,trip.j),blk{trip.nb}.y(trip.i,trip.j));
% probe_legend(length(probe)+1) = string(sprintf('Trip location'));
% 
% 
% legend(probe_plot,probe_legend)
% 
% set(gca,'FontSize',14)
% 
% axis equal

% x = blk.x{6}(:,1,1);
% y = blk.y{6}(:,1,1);
% dx_surf = (x(end)-x(1))/length(x);


%%
figure
hold on
i=3;
skip = 8;
x = blk.x{i};
y = blk.y{i};
x = x(:,:,1);
y = y(:,:,1);
ni = size(x,1);
nj = size(x,2);



for j=[1 skip:skip:ni]
    if (j==1) || (j==ni)
        plot(x(j,:),y(j,:),'r','LineWidth',1)
    else
        plot(x(j,:),y(j,:),'k')
    end
end
for j=[1 skip:skip:nj]
    if (j==1) || (j==nj)
        plot(x(:,j),y(:,j),'b','LineWidth',1)
    else
        plot(x(:,j),y(:,j),'k')
    end
end

%% Cell aspect ratio

figure
hold on
for i=1:NB
    x = blk.x{i};
    y = blk.y{i};
    x = x(:,:,1);
    y = y(:,:,1);

    dxi = x(2:end,:)-x(1:end-1,:);
    dxi(end+1,:) = dxi(end,:);
    dxj = x(:,2:end)-x(:,1:end-1);
    dxj(:,end+1) = dxj(:,end);
    
    dyi = y(2:end,:)-y(1:end-1,:);
    dyi(end+1,:) = dyi(end,:);
    dyj = y(:,2:end)-y(:,1:end-1);
    dyj(:,end+1) = dyj(:,end);
    
    dsi = sqrt(dxi.^2+dyi.^2);
    dsj = sqrt(dxj.^2+dyj.^2);

    min_spacing(i) = min(min(dsi, dsj),[],'all');
    
    aspect_ratio{i} = max(dsi./dsj, dsj./dsi);
    pcolor(blk.x{i},blk.y{i},aspect_ratio{i})
    contour(blk.x{i},blk.y{i},aspect_ratio{i},[11 13 15],'-k')
end
shading('interp')

axis equal
cb = colorbar;

%% O grid analysis

% A = readmatrix('RANS_Ufric_data_Re200k');
% x_RANS = A(:,1);
% [x_RANS, inds, ~] = unique(x_RANS);
% c_RANS = 0.02;
% tau_w = A(inds,2);
% rho = A(inds,3);
% mu = A(inds,4);
% nu = mu./rho;

% Ufric = sqrt(tau_w./rho);

x_grid = [];
y_grid = [];
if NB == 12
    oblocks = [4 6 9 5];
    oblocks_flip = [0, 0, 1, 1];
elseif NB == 9
    oblocks = [3 5 7 4];
    oblocks_flip = [0, 0, 1, 1];
else
    oblocks = 1:length(blk.x);
    oblocks_flip(1:langth(blk.x)) = 0;
end
    
for i=1:length(oblocks)
    x_tmp = blk.x{oblocks(i)}(:,end:-1:1);
    y_tmp = blk.y{oblocks(i)}(:,end:-1:1);
    if oblocks_flip(i) == 1
        x_tmp = x_tmp(end:-1:1,:);
        y_tmp = y_tmp(end:-1:1,:);
    end
    x_grid = [x_grid; x_tmp];
    y_grid = [y_grid; y_tmp];
end
x = x_grid(:,1);
y = y_grid(:,1);
[~,iLE] = min(x);
[~,iTE] = max(x);
iMid = round(0.5*(iLE+iTE));
x = x(iLE:iTE);
tmp = iLE:iTE;
[x, inds, ~] = unique(x);
inds = tmp(inds);
dy = sqrt((x_grid(inds,2)-x_grid(inds,1)).^2 + (y_grid(inds,2) - y_grid(inds,1)).^2);
dx = sqrt((x_grid(inds(2:end),1) - x_grid(inds(1:end-1),1)).^2 + (y_grid(inds(2:end),1) - y_grid(inds(1:end-1),1)).^2);
span = 0.1;
dz = span/(nk-1);

x_cut = x_grid(iMid,:);
y_cut = y_grid(iMid,:);
dy_bl = sqrt((x_cut(2:end)-x_cut(1:end-1)).^2 + (y_cut(2:end) - y_cut(1:end-1)).^2);
expansion_r = dy_bl(2:end)./dy_bl(1:end-1);


% Ufric = interp1(x_RANS,Ufric,x,'pchip','extrap');
% 
% nu = interp1(x_RANS,nu,x,'pchip','extrap');
% nu = nu/c_RANS;
% Xplus = dx.*Ufric(1:end-1)./nu(1:end-1);
% Yplus = dy.*Ufric./nu;
% Zplus = dz.*Ufric./nu;
% 
% figure(3)
% plot(x,Yplus)
% ylabel('First point y^{+}')
% xlabel('x/c')
% xlim([0 1])
% set(gca,'FontSize',16)
% grid on
% 
% figure(4)
% plot(x(1:end-1),Xplus)
% hold on
% plot(x,Zplus)
% xlabel('x/c')
% xlim([0 1])
% legend('\Delta x^{+}','\Delta z^{+}')
% set(gca,'FontSize',16)
% grid on

figure
plot(2:length(expansion_r)+1,expansion_r)
xlabel('y point')
ylabel('O grid expansion ratio')
set(gca,'FontSize',16)

%%


xo = [];
yo = [];
io = [];
jo = [];
blko = [];
for i=1:length(oblocks)
    nb = oblocks(i);
    ni = size(blk.x{nb},1);
    nj = size(blk.x{nb},2);
    xtmp = blk.x{nb}(:,:);
    ytmp = blk.y{nb}(:,:);
    itmp = 1:ni;
    itmp = repmat(itmp',[1 nj]);
    jtmp = 1:nj;
    jtmp = repmat(flip(jtmp), [ni 1]);
    blktmp = nb*ones(ni,nj);
    xtmp = flip(xtmp,2);
    ytmp = flip(ytmp,2);
    if oblocks_flip(i) == 1
        xtmp = flip(xtmp);
        ytmp = flip(ytmp);
        itmp = flip(itmp);
    end
    if size(xo,1) == 0
        xo = xtmp; yo = ytmp;
        io = itmp; jo = jtmp;
        blko = blktmp;
    else
        xo = [xo; xtmp(2:end,:)];
        yo = [yo; ytmp(2:end,:)];
        io = [io; itmp(2:end,:)];
        jo = [jo; jtmp(2:end,:)];
        blko = [blko; blktmp(2:end,:)];
    end
    
    if xo(1,1) == xo(end,1)
        xo = xo(1:end-1,:);
        yo = yo(1:end-1,:);
        io = io(1:end-1,:);
        jo = jo(1:end-1,:);
        blko = blko(1:end-1,:);
    end

end
xsurf = xo(:,1);
[~, iLE] = min(xsurf);
[~, iTE] = max(xsurf);

ssurf = zeros(1,iTE-iLE+1);
for i = iLE:iTE
    if i>iLE
        dx = xo(i,1) - xo(i-1,1);
        dy = yo(i,1) - yo(i-1,1);
        ds = sqrt(dx^2 + dy^2);
        ssurf(i+1-iLE) = ssurf(i-iLE) + ds;
    end
end
c = sqrt((xo(iLE,1)-xo(iTE,1))^2 + (yo(iLE,1)-yo(iTE,1))^2);
yp1 = Re2offset(Re,1,c);

dx = xo(2:end,1) - xo(1:end-1,1);
dy = yo(2:end,1) - yo(1:end-1,1);
ds = sqrt(dx.^2 + dy.^2);
ds = ds(iLE:iTE);

dx = xo(:,2) - xo(:,1);
dy = yo(:,2) - yo(:,1);
dn = sqrt(dx.^2 + dy.^2);
dn = dn(iLE:iTE);
xnow = xo(iLE:iTE);

figure
yyaxis left
plot(xnow, dn/yp1);
yyaxis right
plot(xnow, ds/yp1);
hold on
yline(dz/yp1, '--')

%%

fprintf('%d points in i-j plane\n',nij_pts)
fprintf('%d processors\n', nk*sum(prod(procdims,2))/npp)
fprintf('Total points: %d\n', nij_pts*nk);
fprintf('Min spacing: %6.4e\n', min(min_spacing));
end