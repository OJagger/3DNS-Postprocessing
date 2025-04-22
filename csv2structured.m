function s = csv2structured(basecase, fname, metadata2set, sample_blk, iPS)
    % tecfile = '/Data/ojj23/dns-work/channel-flows/RANS/M1-2_Re2k2-bsl/tecplot.dat';
    %[zone, VARlist] = tec2mat(tecfile,'debug');

    if nargin < 5
        iPS = false;
    end

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
    

    s = RANSSlice(blk,gas,bcs,[],[],iPS);
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
    names = string(data.Properties.VariableNames);
    if ismember("x_coordinate",names)
        names = struct('x','x_coordinate', ...
            'y','y_coordinate', ...
            'ro','density', ...
            'u','x_velocity', ...
            'v','y_velocity', ...
            'p', 'pressure', ...
            'k', 'turb_kinetic_energy', ...
            'om', 'specific_diss_rate', ...
            'sp', 'modified_viscosity', ...
            'mut', 'viscosity_turb', ...
            'st', 'strain_rate_mag', ...
            'walldist', ["walldist", "cell_wall_dist", "cell_wall_distance"]);
  
    else
        names = struct('x','Points_0', ...
            'y','Points_1', ...
            'ro','density', ...
            'u','relativeVelocityVector_0', ...
            'v','relativeVelocityVector_1', ...
            'p', 'staticPressure', ...
            'k', 'turbulentKineticEnergy', ...
            'om', 'turbulentOmega', ...
            'sp', 'spalartVariable', ...
            'mut', 'turbulentViscosity', ...
            'st', 'strain_rate_mag', ...
            'walldist', ["wallDistance"]);
    end

    data_fields = convertCharsToStrings(data.Properties.VariableNames);

    xv = data.(names.x);
    yv = data.(names.y);
    rv = data.(names.ro);
    uv = data.(names.u);
    vv = data.(names.v);
    wv = zeros(size(uv));
    pv = data.(names.p);

    ri = scatteredInterpolant(xv, yv, rv,'linear','boundary');
    ui = scatteredInterpolant(xv, yv, uv,'linear','boundary');
    vi = scatteredInterpolant(xv, yv, vv,'linear','boundary');
    wi = scatteredInterpolant(xv, yv, wv,'linear','boundary');
    pi = scatteredInterpolant(xv, yv, pv,'linear','boundary');

    save_k = false;
    if ismember(string(names.k), data_fields)
        ki = scatteredInterpolant(xv,yv,data.(names.k),'linear','boundary');
        save_k = true;
    end

    save_om = false;
    if ismember(string(names.om), data_fields)
        oi = scatteredInterpolant(xv,yv,data.(names.om),'linear','boundary');
        save_om = true;
    end

    save_mut = false;
    if ismember(string(names.mut), data_fields)
        mi = scatteredInterpolant(xv,yv,data.(names.mut),'linear','boundary');
        save_mut = true;
    end

    save_st = false;
    if ismember(string(names.st), data_fields)
        si = scatteredInterpolant(xv,yv,data.(names.st),'linear','boundary');
        save_st = true;
    end

    save_walldist = false;
    iname = find(ismember(string(names.walldist), data_fields));
    if length(iname) > 0
        di = scatteredInterpolant(xv,yv,data.(names.walldist(iname(1))),'linear','boundary');
        save_walldist = true;
    end
    % elseif ismember(string(names.walldist2), data_fields)
    %     di = scatteredInterpolant(xv,yv,data.(names.walldist2),'linear','boundary');
    %     save_walldist = true;
    % elseif ismember(string(names.walldist3), data_fields)
    %     di = scatteredInterpolant(xv,yv,data.(names.walldist3),'linear','boundary');
    %     save_walldist = true;
    % end

    save_spalart = false;
    if ismember(string(names.sp), data_fields)
        spi = scatteredInterpolant(xv,yv,data.(names.sp),'linear','boundary');
        save_spalart = true;
    end

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
            s.turb_model = 'ko';
            disp('Interpolating omega')
            s.omega{ib} = oi(blk.x{ib}, blk.y{ib});
        end

        if save_mut
            disp('Interpolating mut')
            s.mut_store{ib} = mi(blk.x{ib}, blk.y{ib});
        end

        if save_st
            disp('Interpolating S')
            s.StR_store{ib} = si(blk.x{ib}, blk.y{ib});
        end

        if save_walldist
            disp('Interpolating wall distance')
            s.blk.walldist{ib} = di(blk.x{ib}, blk.y{ib});
        end

        if save_spalart
            s.turb_model = 'sa';
            disp('Interpolating Spalart variable');
            s.sp{ib} = spi(blk.x{ib}, blk.y{ib});
        end


    end
    i = 1;
    while sum(isnan(s.u{1}(i,:))) > 0
        i = i+1;
    end
    fprintf('Inlet index for inlet stagnation conditions sampling: %d\n', i)
    s.getBCs([],(i+40:i+100));
end
    