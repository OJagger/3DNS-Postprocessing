function sol = read_mrchbl_output(path, pitch)

    if nargin<2
        scale = 1;
    else
        f = fileparts(path);
        f = fopen(fullfile(f, 'blade.mises'));
        l = fgetl(f);
        temp = str2num(char(split(fgetl(f))));
        scale = pitch/temp(5);
    end


    system([char("sed -i 's/\*\*\*\*\*\*\*\*\*\*/       NaN/g' ") path])
    f = fopen(path, 'r');

    for i=1:3
        fgetl(f);
    end
    temp = split(fgetl(f), '|');
    temp = str2num(char(split(strtrim(temp(1)))));
    sol.Klag = temp(1);
    sol.Kp = temp(2);
    sol.Kd = temp(3);
    l = fgetl(f);
    if contains(l, '|')
        temp = split(l, '|');
        temp = str2num(char(split(strtrim(temp(1)))));
        sol.xShock = temp(1);
        sol.fctshk = temp(2);
        fgetl(f);
        fgetl(f);
    else
        fgetl(f);
    end

    A = fscanf(f, '%f', [18 Inf]);

    fclose(f);

    toin = 300;
    gam = 1.4;
    cp = 1005;
    rgas = cp*(gam-1)/gam;
    a0 = sqrt(gam*rgas*toin);

    sol.x = A(1,:)*scale;
    sol.Uea0 = A(2,:);
    sol.Ue = sol.Uea0*a0;
    sol.M = M_VT0(sol.Ue, toin, gam, cp);
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