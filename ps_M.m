function ps = ps_Mp0(M, p0, gam)
    
    if nargin < 3
        gam = 1.4;
    end

    ps = p0*p_p0(M, gam)*(1+2*gam*(M^2-1)/(gam+1));

end