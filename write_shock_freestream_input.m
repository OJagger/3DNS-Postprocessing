function write_shock_freestream_input(dir, Min, xshock, Lshock)

gam = 1.4;
cp = 1005;
rgas = cp*(gam-1)/gam;
Poin = 1e5;
Toin = 300;

fM = 1+0.5*(gam-1)*Min^2;
pin = Poin*fM^(-gam/(gam-1));
tin = Toin/fM;
roin = pin/(rgas*tin);
vin = Min*sqrt(gam*rgas*tin);
Etin = pin/(gam-1) + 0.5*roin*vin^2;

% Post shock conditions
Ms = sqrt(fM/(gam*Min^2 - 0.5*(gam-1)));
ps = pin*(1+2*gam*(Min^2-1)/(gam+1));
ros = 0.5*roin*(gam+1)*Min^2/fM;
Ts = ps/(ros*rgas);
vs = Ms*sqrt(gam*rgas*Ts);
Ets = ps/(gam-1) + 0.5*ros*vs^2;

x = [linspace(0,(xshock-Lshock/2), 4) xshock linspace(xshock+Lshock/2,1, 4)];
p = zeros(size(x));
p(1:4) = pin;
p(end-4:end)=ps;
p(5) = 0.5*(pin+ps);

v = zeros(size(x));
v(1:4) = vin;
v(end-4:end)=vs;
v(5) = 0.5*(vin+vs);

path = fullfile(dir, 'freestream.txt');
f = fopen(path, 'w');
fprintf(f, '%d\n', length(x));
for i=1:length(x)
    fprintf(f, '%f %f %f\n', [x(i), v(i), p(i)]);
end
fclose(f)
fprintf('vin: %f\n', vin)
fprintf('pexit: %f\n', ps)
end