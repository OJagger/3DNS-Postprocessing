function s = HYDRA_csv2structured(basecase, fname, beta1_fac)
    tecfile = '/Data/ojj23/dns-work/channel-flows/RANS/M1-2_Re2k2-bsl/tecplot.dat';
    [zone, VARlist] = tec2mat(tecfile,'debug');

    if nargin < 3
        beta1 = [];
    end

    blk = basecase.blk;
    gas = basecase.gas;
    bcs = basecase.bcs;
    

    s = RANSSlice(blk,gas,bcs);
    s.blk = basecase.blk;
    s.gas = basecase.gas;
    s.beta1_fac = beta1_fac;
    gam  = gas.gam;

    data = readtable(fname);

    xv = data.Points_0;
    yv = data.Points_1;
    rv = data.density;
    uv = data.relativeVelocityVector_0;
    vv = data.relativeVelocityVector_1;
    wv = data.relativeVelocityVector_2;
    pv = data.staticPressure;
    kv = data.turbulentKineticEnergy;
    ov = data.turbulentOmega;
    dv = data.wallDistance;

    ri = scatteredInterpolant(xv, yv, rv);
    ui = scatteredInterpolant(xv, yv, uv);
    vi = scatteredInterpolant(xv, yv, vv);
    wi = scatteredInterpolant(xv, yv, wv);
    pi = scatteredInterpolant(xv, yv, pv);
    ki = scatteredInterpolant(xv, yv, kv);
    oi = scatteredInterpolant(xv, yv, ov);
    di = scatteredInterpolant(xv, yv, dv);

    for ib = 1:length(blk.x)
        fprintf('Interpolating block %d/%d\n',[ib, length(blk.x)]);
        disp('Interpolating ro')
        s.ro{ib} = ri(blk.x{ib}, blk.y{ib});
        disp('Interpolating u')
        s.u{ib} = ui(blk.x{ib}, blk.y{ib});
        disp('Interpolating v')
        s.v{ib} = vi(blk.x{ib}, blk.y{ib});
        disp('Interpolating w')
        s.w{ib} = wi(blk.x{ib}, blk.y{ib});
        disp('Interpolating p')
        pnow = pi(blk.x{ib}, blk.y{ib});
        s.Et{ib} = pnow/(gam-1) + 0.5*s.ro{ib}.*(s.u{ib}.^2 + s.v{ib}.^2 + s.w{ib}.^2);
        disp('Interpolating k')
        s.k{ib} = ki(blk.x{ib}, blk.y{ib});
        disp('Interpolating omega')
        s.omega{ib} = oi(blk.x{ib}, blk.y{ib});
        disp('Interpolating wallDistance')
        s.blk.walldist{ib} = di(blk.x{ib}, blk.y{ib});

    end
    s.getBCs;
end
    