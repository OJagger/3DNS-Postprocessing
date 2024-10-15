function [ydash] = MISES_ydash(x, y, xsurf, Me, Ue, nue, data)

    Me = reshape(Me, [], 1);
    % Ue = reshape(Ue, [], 1);

    ydash = zeros(3,1);

    dUdx = (Ue(3:end)-Ue(1:end-2))./(xsurf(3:end)-xsurf(1:end-2));
    dUdx = [(Ue(2)-Ue(1))./(xsurf(2)-xsurf(1)); dUdx; ...
        (Ue(end)-Ue(end-1))./(xsurf(end)-xsurf(end-1))];
    dUdx = interp1(xsurf, dUdx, x);

    dMdx = (Me(3:end)-Me(1:end-2))./(xsurf(3:end)-xsurf(1:end-2));
    dMdx = [(Me(2)-Me(1))./(xsurf(2)-xsurf(1)); dMdx; ...
        (Me(end)-Me(end-1))./(xsurf(end)-xsurf(end-1))];
    dMdx = interp1(xsurf, dMdx, x);

    Uea = interp1(xsurf, Ue, x);
    nuea = interp1(xsurf, nue, x);
    Ma = interp1(xsurf, Me, x);

    if isfield(data, "theta")
        theta = interp1(xsurf, data.theta, x);
    else
        theta = y(1);
    end

    Reth = theta*Uea/nuea;

    if isfield(data, "Hs")
        Hs = interp1(xsurf, data.Hs, x);
    else
        Hs = y(2);
    end

    Hks = MISES_correlations.fHks_Hs(Hs, Ma);

    if isfield(data, "Hk")
        Hk = interp1(xsurf, data.Hk, x);
    else
        Hk = MISES_correlations.fHk_Hks(Hks, Reth);
    end

    if isfield(data, "H")
        H = interp1(xsurf, data.H, x);
    else
        H = MISES_correlations.fH(Ma, Hk);
    end

    if isfield(data, "delStar")
        delStar = interp1(xsurf, data.delStar, x);
        H = delStar/theta;
    else
        delStar = H*theta;
    end

    % Hks = MISES_correlations.fHks(Hk, Reth);

    % if isfield(data, "Hs")
    %     Hs = interp1(xsurf, data.Hk, x);
    % else
    %     Hs = MISES_correlations.fHs(Ma, Hks);
    % end

    Hss = MISES_correlations.fHss(Ma, Hk);

    if isfield(data, "cf")
        cf = interp1(xsurf, data.cf, x);
    else
        cf = MISES_correlations.fCf(Hk, Reth, Ma);
    end
    
    Us = MISES_correlations.fUs(Hs,Hk,H);
    Cd = MISES_correlations.fCd(cf, Us, y(3));
    CtEQ = MISES_correlations.fCtEQ(Hs, Us, H, Hk);
    del = MISES_correlations.fDel(theta, Hk, delStar);

    ydash(1) = cf/2 - (H + 2 - Ma^2)*(theta/Uea)*dUdx;

    % dHsdx = 2*Cd/theta - Hs*cf/(2*theta) - dUdx*(2*Hss + Hs*(1-H))/Uea;
    % 
    % delThick = 5e-12;                                                           % Change in quantities for differencing
    % delM = 1e-3;
    % 
    % Hs1 = MISES_correlations.fHs_ds(delStar+delThick, theta, Ma, Reth);
    % Hs0 = MISES_correlations.fHs_ds(delStar-delThick, theta, Ma, Reth);
    % dHs_dds = (Hs1-Hs0)/(2*delThick);                                          % Partial derivative of Hs wrt delta star.
    % 
    % Hs1 = MISES_correlations.fHs_ds(delStar, theta+delThick, Ma, Reth);
    % Hs0 = MISES_correlations.fHs_ds(delStar, theta-delThick, Ma, Reth);
    % dHs_dtheta = (Hs1-Hs0)/(2*delThick);                                       % Partial derivative of Hs wrt theta.
    % 
    % Hs1 = MISES_correlations.fHs_ds(delStar, theta, Ma+delM, Reth);
    % Hs0 = MISES_correlations.fHs_ds(delStar, theta, Ma-delM, Reth);
    % dHs_dMe = (Hs1-Hs0)/(2*delM);                                          % Partial derivative of Hs wrt Me.
    % 
    % ydash(2) = (dHsdx - dHs_dtheta*ydash(1) - dHs_dMe*dMdx)/dHs_dds;

    ydash(2) = (2*Cd - Hs*cf/2 - (2*Hss+Hs*(1-H))*(theta/Uea)*dUdx) ...
        / theta;

    % Set Kcorr value to default or custom
    if isfield(data, "Kcorr")
        Kcorr = data.Kcorr;
    else
        Kcorr = 4.2;
    end

    % Set Kd value if included
    if isfield(data, "Kd")
        Kd = data.Kd;
    else
        Kd = 0;
    end

    % Set Kp value if included
    if isfield(data, "Kp")
        Kp = data.Kp;
    else
        Kp = 0;
    end

    A = 6.7; B = 0.75; % MISES Eqm locus consts


    ydash(3) = (y(3)/del)*( ...
        Kcorr*(sqrt(CtEQ) - sqrt(y(3))) ...
        + Kd*2*del*H/(B*theta)*(cf/2 - ((Hk - 1)/(A*Hk))^2) ...
        - Kp*2*(del/Uea)*dUdx);

    if x>0.3716
        disp('');
    end
end