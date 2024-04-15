function val = ro_ro0(M, gam)

if nargin<2
    gam=1.4;
end

val = (1 + 0.5*(gam-1)*M.^2).^(-1/(gam-1));

end