function sol = read_mrchbl_output(path)

    f = fopen(path, 'r');

    for i=1:6
        fgetl(f);
    end

    A = fscanf(f, '%f', [18 Inf]);

    fclose(f);

    sol.x = A(1,:);
    sol.Ue = A(2,:);
    sol.delStar = A(5,:);
    sol.theta = A(6,:);
    sol.H = A(7,:);
    sol.H_k = A(8,:);
    sol.H_bar = A(9,:);
    sol.H_ke = A(10,:);
    sol.H_rho = A(11,:);
    sol.Cf = A(12,:);
    sol.cd = A(13,:);
    sol.ct = A(14,:);
    sol.Us = A(15,:);
    sol.Uq = A(16,:);
    sol.cteq = A(17,:);
    sol.Re_theta = A(18,:);
    sol.Pr = (1-sol.Us).*sol.ct;
    sol.Pr_eq = (1-sol.Us).*sol.cteq;
    

end