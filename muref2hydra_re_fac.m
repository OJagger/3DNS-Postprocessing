function refac = muref2hydra_re_fac(gas)

    lref = 1;
    roref = 1.226;
    pref = 1.013e5;
    Tref = 288;
    uref = sqrt(pref/roref);
    muref = 5.0872e-8;

    mu0 = muref*roref*uref*lref;
    mu = sutherland_mu(Tref, gas.mu_ref, gas.mu_cref, gas.mu_tref);

    refac = mu/mu0;

end
