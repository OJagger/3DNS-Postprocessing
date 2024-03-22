function s = Fluent_csv2structured(basecase, fname, metadata2set, sample_blk)
    % tecfile = '/Data/ojj23/dns-work/channel-flows/RANS/M1-2_Re2k2-bsl/tecplot.dat';
    %[zone, VARlist] = tec2mat(tecfile,'debug');


    if nargin < 4 || isempty(sample_blk)
        blk = basecase.blk;
    else
        blk = sample_blk;
        for var = ["inlet_blocks", "viewarea", "oblocks", "oblocks_flip", "aspect"]
            blk.(var) = basecase.blk.(var);
        end
    end

    gas = basecase.gas;
    bcs = basecase.bcs;
    

    s = RANSSlice(blk,gas,bcs);
    if nargin > 3 && ~isempty(metadata2set)
        keys = fieldnames(metadata2set);
        for i=1:length(keys)
            s.(keys(i)) = metadata2set.(keys(i));
        end
    end

%     s.blk = basecase.blk;
%     s.gas = basecase.gas;
%     s.beta1_fac = beta1_fac;
    gam  = gas.gam;

    data = readtable(fname);

    xv = data.x_coordinate;
    yv = data.y_coordinate;
    rv = data.density;
    uv = data.x_velocity;
    vv = data.y_velocity;
    wv = zeros(size(uv));
    pv = data.pressure;
%     kv = data.turb_kinetic_energy;
%     ov = data.specific_diss_rate;
%     dv = data.wallDistance;

    ri = scatteredInterpolant(xv, yv, rv,'linear','none');
    ui = scatteredInterpolant(xv, yv, uv,'linear','none');
    vi = scatteredInterpolant(xv, yv, vv,'linear','none');
    wi = scatteredInterpolant(xv, yv, wv,'linear','none');
    pi = scatteredInterpolant(xv, yv, pv,'linear','none');
%     ki = scatteredInterpolant(xv, yv, kv,'linear','none');
%     oi = scatteredInterpolant(xv, yv, ov,'linear','none');

    save_k = false;
    if ismember("turb-kinetic-energy", convertCharsToStrings(data.Properties.VariableNames))
        ki = scatteredInterpolant(xv,yv,data.turb-kinetic-energy,'linear','none');
        save_k = true;
    end

    save_om = false;
    if ismember("specific_diss_rate", convertCharsToStrings(data.Properties.VariableNames))
        oi = scatteredInterpolant(xv,yv,data.specific_diss_rate,'linear','none');
        save_om = true;
    end

    save_mut = false;
    if ismember("viscosity_turb", convertCharsToStrings(data.Properties.VariableNames))
        mi = scatteredInterpolant(xv,yv,data.viscosity_turb,'linear','none');
        save_mut = true;
    end

    save_mut = false;
    if ismember("viscosity_turb", convertCharsToStrings(data.Properties.VariableNames))
        mi = scatteredInterpolant(xv,yv,data.viscosity_turb,'linear','none');
        save_mut = true;
    end
%     save_pr = false;
%     if ismember("production_of_k", convertCharsToStrings(data.Properties.VariableNames))
%         pri = scatteredInterpolant(xv,yv,data.production_of_k,'linear','none');
%         save_pr = true;
%     end
    save_st = false;
    if ismember("strain_rate_mag", convertCharsToStrings(data.Properties.VariableNames))
        si = scatteredInterpolant(xv,yv,data.strain_rate_mag,'linear','none');
        save_st = true;
    end
%     di = scatteredInterpolant(xv, yv, dv);

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
        if save_k
            disp('Interpolating k')
            s.k{ib} = ki(blk.x{ib}, blk.y{ib});
        end
        if save_om
            disp('Interpolating omega')
            s.omega{ib} = oi(blk.x{ib}, blk.y{ib});
        end
        if save_mut
            disp('Interpolating mut')
            s.mut_store{ib} = mi(blk.x{ib}, blk.y{ib});
        end
%         if save_pr
%             disp('Interpolating pr')
%             s.Pr_store{ib} = pri(blk.x{ib}, blk.y{ib});
%         end
        if save_st
            disp('Interpolating S')
            s.StR_store{ib} = si(blk.x{ib}, blk.y{ib});
        end
%         disp('Interpolating wallDistance')
%         s.blk.walldist{ib} = di(blk.x{ib}, blk.y{ib});

    end
    i = 1;
    while sum(isnan(s.u{1}(i,:))) > 0
        i = i+1;
    end
    fprintf('Inlet index for inlet stagnation conditions sampling: %d\n', i)
    s.getBCs([],(i+40:i+100));
end
    