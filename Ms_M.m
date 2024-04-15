function Ms = Ms_M(M, gam)
    
    if nargin < 2
        gam = 1.4;
    end

    fM = 1+0.5*(gam-1)*M^2;
    Ms = sqrt(fM/(gam*M^2 - 0.5*(gam-1)));

end