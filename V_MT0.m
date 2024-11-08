function V = V_MT0(M,T0,gam,cp)
    
    if nargin < 3 || isempty(gam)
        gam = 1.4;
    end
    if nargin < 4 || isempty(cp)
        cp = 1005;
    end

    T = T0*T_T0(M,gam);
    V = sqrt( 2*cp * (T0-T) );
    
end
