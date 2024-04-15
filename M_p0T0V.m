function M = M_p0T0V(p0,T0,v,gam,cp)

    if nargin < 4 || isempty(gam)
        gam = 1.4;
    end
    if nargin < 5 || isempty(cp)
        cp = 1005;
    end

    rgas = cp*(gam-1)/gam;

    T = T0 - v.^2/(2*cp);
    M = v./sqrt(gam*rgas*T);
end