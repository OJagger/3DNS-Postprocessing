function newcase = setup_channel_case(name, path, blk, Min, theta, varargin)

    p = inputParser;

    addRequired(p,'name');
    addRequired(p,'path');
    addRequired(p,'blk')
    addRequired(p,'Min');
    addRequired(p,'theta');
    addParameter(p,'Reth',[]);
    addParameter(p,'mu_ref',[])
    addParameter(p,'xshock',[]);

    
    parse(p, name, path, blk, Min, theta, varargin{:})

    bcs.Poin = 1e5;
    bcs.Toin = 300;
    bcs.alpha = 0;
    bcs.gamma = 0;
    bcs.aturb = 150;
    bcs.nturb = 100;
    bcs.iradprof = 0;
    bcs.g_z = 0;
    bcs.theta = theta;

    solver.nkprocs = 1;
    solver.nk = blk.nk;
    solver.niter = 1000;
    solver.nwrite = 10000;
    solver.ncut = 1000;
    solver.cfl = 1;
    solver.sigma = 0.0300;
    solver.span = 0.1000;
    solver.fexpan = 1;
    solver.irestart = 0;
    solver.istats = 0;
    solver.istability = 0;
    solver.ilam = 1;
    solver.npp = blk.nk;


    gas.gam = 1.4;
    gas.cp = 1005;
    gas.mu_tref = 273;
    gas.mu_cref = 110.4000;
    gas.pr = 0.7200;
    rgas = gas.cp*(gas.gam-1)/gas.gam;

    fM = 1+0.5*(gas.gam-1)*Min^2;
    pin = bcs.Poin*fM^(-gas.gam/(gas.gam-1));
    Tin = bcs.Toin*T_T0(Min, gas.gam);
    roin = pin/(rgas*Tin);
    bcs.vin = sqrt(2*gas.cp*(bcs.Toin-Tin));
    Ms = sqrt(fM/(gas.gam*Min^2 - 0.5*(gas.gam-1)));
    ps = pin*(1+2*gas.gam*(Min^2-1)/(gas.gam+1));

    bcs.pexit = ps;

    if isempty(p.Results.Reth)
        if isempty(p.Results.mu_ref)
            fprintf('Must specify either Reth mu_ref at inlet\n');
            return
        else
            gas.mu_ref = p.Results.mu_ref;         
        end
    else
        mu_in = theta*bcs.vin/p.Results.Reth;
        gas.mu_ref = sutherland_mu_ref(mu_in, Tin, gas.mu_cref, gas.mu_tref);
    end

    newcase = DNS_channel;
    newcase.casename = name;
    newcase.blk = blk;
    newcase.gas = gas;
    newcase.solver = solver;
    newcase.bcs = bcs;
    newcase.casepath = fullfile(path, name);
    newcase.casetype = 'gpu';

end