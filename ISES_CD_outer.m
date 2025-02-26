function [Pr, Ct, CtEQ, pgterm, diffterm, Us] = ISES_CD_outer(s, theta, H_ke, H, H_k, delStar, Ue, p_wall, ro, cf, istart, Kcorr, Kp, Kd)

    s = reshape(s,[],1);
    theta = reshape(theta,[],1);
    H_ke = reshape(H_ke,[],1);
    H = reshape(H,[],1);
    H_k = reshape(H_k,[],1);
    delStar = reshape(delStar,[],1);
    Ue = reshape(Ue,[],1);
    cf = reshape(cf,[],1);
    p_wall = reshape(p_wall,[],1);
    ro = reshape(ro,[],1);

    if nargin < 10 || isempty(Kcorr)
        Kcorr= 4.2;
    end

    Kcorr;

    if nargin < 11 || isempty(Kp)
        Kp = 0;
    end

    if nargin < 12 || isempty(Kd)
        Kd = 0;
    end

    A = 6.7; B = 0.75; % MISES Eqm locus consts
    
    
    dUedx = (Ue(3:end) - Ue(1:end-2))./(s(3:end) - s(1:end-2));
    dUedx = [dUedx(1); dUedx; dUedx(end)];

    dpdx = (p_wall(3:end) - p_wall(1:end-2))./(s(3:end) - s(1:end-2));
    dpdx = [dpdx(1); dpdx; dpdx(end)];


    Us = 0.5*H_ke.*(1-4*(H_k-1)./(3*H));
    del = theta.*(3.15+1.72./(H_k - 1)) + delStar;
    CtEQ = H_ke.*(0.015./(1-Us)).*(H_k-1).^3./(H_k.^2.*H);
   % CtEQ = 0.02456*((H_k-1)./H_k).^3./(1-Us);
    Ct(1:istart-1) = NaN;
    Ct(istart) = CtEQ(istart);


    for i=istart:length(s)-1
        dCt_ds = Kcorr*(sqrt(CtEQ(i)) - sqrt(Ct(i)));

        diffterm(i) = (2*del(i)*H(i)/(B*theta(i)))*(cf(i)/2 - ((H_k(i) - 1)/(A*H_k(i)))^2);
        dCt_ds = dCt_ds + Kd*diffterm(i);

%         pgterm(i) = (2*del(i)/Ue(i))*dUedx(i);
        pgterm(i) = (2*del(i)/Ue(i)^2/ro(i))*dpdx(i);
        dCt_ds = dCt_ds + Kp*pgterm(i);
        
        dCt_ds = dCt_ds*Ct(i)/del(i);
        Ct(i+1) = Ct(i) + dCt_ds*(s(i+1)-s(i));
    end

    diffterm(end+1) = NaN;
    pgterm(end+1) = NaN;


    diffterm = reshape(diffterm,[],1);
    pgterm = reshape(pgterm, [], 1);
    Ct = reshape(Ct,[],1);
    Us = reshape(Us,[],1);

    Pr = Ct.*(1-Us);


end