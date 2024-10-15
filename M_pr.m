function M = M_pr(pr, gam)

if nargin < 2
    gam = 1.4;
end

M = sqrt(2*(pr.^(-(gam-1)/gam)-1)/(gam-1));

end