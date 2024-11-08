function sol = read_mrchbl_output(path)

    f = fopen(path, 'r');

    for i=1:6
        fgetl(f);
    end

    A = fscanf(f, '%f', [15 Inf]);

    fclose(f);

    sol.x = A(1,:);
    sol.Ue = A(2,:);
    sol.delStar = A(5,:);
    sol.theta = A(6,:);
    sol.H = A(7,:);
    sol.H_k = A(8,:);
    sol.H_bar = A(9,:);
    sol.H_ke = A(10,:);
    sol.Cf = A(11,:);
    sol.cd = A(12,:);
    sol.ct = A(13,:);
    sol.Us = A(14,:);
    sol.Re_theta = A(15,:);
    sol.Pr = (1-sol.Us).*sol.ct;
    

end