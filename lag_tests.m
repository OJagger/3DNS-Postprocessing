clear

% boundary solver
%rey = 3e5;
rey = 2e6;

% transition point Re_theta;
%retheta_trans = 200;
retheta_trans = 400;

% Lag constant
Klag = 2.8;%5.6;

% streamwise points
ni = 20000;
x = logspace(-4,0,ni);
 
ue(1:ni) = 1.0;

step = 0; %0.2*(1+tanh((x-0.5)*500))*0.5;
ue = ue - step;

% initial conditions
h = 2.59;
hstar = 1.573;

dstar(1) = 1.72d0*sqrt(x(1))/sqrt(rey);     
theta(1) = dstar(1)/h;
eps(1) = theta(1)*hstar;

ue2 = ue.*ue;
ue3 = ue2.*ue;

logue = log(ue);

itrans = 0;

for i=2:ni
ueav = (ue(i) + ue(i-1))*0.5;
ue2av = (ue2(i) + ue2(i-1))*0.5;
ue3av = (ue3(i) + ue3(i-1))*0.5;
   
rt = rey.*theta(i-1).*ue(i-1);
rex = rey*x(i)*ue(i);
h_old = h;

if(rt<retheta_trans)% laminar

% solve for h    
for iter=1:100
hstar_guess = 1.515 + 0.076*((4-h).^2)./h;
dh = (hstar - hstar_guess)/(-0.076*( (4-h).^2 + 2*h*(4-h))/(h*h));
if(abs(dh)<1e-4)
    break
end
h = h + dh;
end

cf = (-0.067 + 0.01977*((7.4 - h).^2)./(h-1)).*(2./rt) ;   
cd = (0.207 + 0.00205*((4 - h).^5.5)).*hstar./(2*rt);
ctau = 0.0;
ctaueq = 0.0;
us = 0.5*hstar.*(1 - (4/3)*(h-1)./h);

dlogue_eq = 0;
del = 5.0*x(i)/sqrt(rex);

else % turbulent
    
del = theta(i-1).*(3.15 + 1.72./(h-1)) + dstar(i-1);    

ho = 4;
if(rt>400)
    ho = 3 + (400/rt);
end

% solve for h
for iter=1:100
if(h<ho)   
hstar_guess = 1.505 + (4)./rt + (0.165 - (1.6)./sqrt(rt))*((ho-h).^1.6)./h;
else
hstar_guess = 1.505 + (4./rt) + ((h-ho).^2).*( (0.04./h) + 0.007*log(rt)./((h - ho + (4./log(rt))).^2) );    
end
dh = (hstar - hstar_guess)/(-0.076*( (4-h).^2 + 2*h*(4-h))/(h*h));
if(abs(dh)<1e-4)
    break
end
h = h + dh;
end

us = 0.5*hstar.*(1 - (4/3)*(h-1)./h);
ctaueq = (0.015*hstar./(1-us)).*( (h-1).^3 )./(h.*(h.^2)) ;    
%ctau = ctaueq;

if(itrans == 0)
ctau = (3.24*ctaueq.*exp(-6.6./(h-1)));
itrans = 1;
dctau =  ctau*( (sqrt(ctaueq) - sqrt(ctau) )*Klag/del ...
+ (8/(3*dstar(i-1)))*(cf*0.5 - ((h-1)/(6.7*h))^2) ...    
- 2*dlogue );
end

% shear-lag model
dctau_old = dctau;
dctau =  ctau*( (sqrt(ctaueq) - sqrt(ctau) )*Klag/del ...
+ (8/(3*dstar(i-1)))*(cf*0.5 - ((h-1)/(6.7*h))^2) ...    
- 2*dlogue );

dlogue_eq = -(1/2)*(8/(3*dstar(i-1)))*(cf*0.5 - ((h-1)/(6.7*h))^2);

ctau = ctau + (dctau*1.5 - 0.5*dctau_old)*dx;

cf = 0.3*exp(-1.33*h).*(log10(rt).^(-1.74 - 0.31*h)) ...
    + 0.00011*(tanh(4 - (h/0.875)) - 1);
cd = cf*0.5*us + ctau*(1-us);

end

dx = x(i)-x(i-1);    
due2 = (ue2(i)-ue2(i-1))/dx;
due3 = (ue3(i)-ue3(i-1))/dx;
dlogue = (logue(i)-logue(i-1))/dx;
dthet = cf*0.5 - (2 + h)*theta(i-1)*dlogue;

dhstar = (2*cd - hstar*cf*0.5 - hstar*(1-h)*theta(i-1)*dlogue)./theta(i-1);

if i>2
theta(i) = theta(i-1) + (dthet*1.5 - dthet_old*0.5)*dx;
hstar = hstar + (dhstar*1.5 - dhstar_old*0.5)*dx;
else
theta(i) = theta(i-1) + dthet*dx;
hstar = hstar + dhstar*dx;
end

dthet_old = dthet;
cf_old = cf;
cd_old = cd;
hstar_old = hstar;
dhstar_old = dhstar;

eps(i) = theta(i)*hstar;
dstar(i) = theta(i)*h;


Reth(i) = rt;
Re_x(i) = rex;
CF(i) = cf;
H(i) = h;
Pr(i) = ctau*(1-us);
Pr_eq(i) = ctaueq*(1-us);
Terms.Dlogue_eq(i) = -dlogue_eq*del;
Terms.Dlogue(i) = -dlogue*del;
Terms.Ctau(i) = sqrt(ctau)*Klag;
Terms.Ctau_eq(i) = sqrt(ctaueq)*Klag;

end

figure(5)
plot(H,Pr,'k.')
hold on
plot(H,Pr_eq,'-r')
ylabel('Pr'),xlabel('H')

figure(6)
plot(x,Pr,'-k.',x,Pr_eq,'r')
ylabel('Pr'),xlabel('x')

figure(7)
plot(x,Terms.Dlogue,'-k.',x,Terms.Dlogue_eq,'r')
ylabel('Pressure gradient terms'),xlabel('x')

figure(8)
plot(x,Terms.Ctau,'-k.',x,Terms.Ctau_eq,'r')
ylabel('Lag terms'),xlabel('x')
