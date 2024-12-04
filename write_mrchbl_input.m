function write_mrchbl_input(fpath, x, Ue_tilde, M, Reth, xtrip, th, ds, ct, Kcorr, Kp, Kd, xshock, fctshk)
    
    if nargin < 10 || isempty(Kcorr)
        Kcorr = 5.6;
    end
    if nargin < 11 || isempty(Kp)
        Kp = 1.0;
    end
    if nargin < 12 || isempty(Kd)
        Kd = 1.0;
    end
    if nargin < 14 || isempty(xshock) || isempty(fctshk)
        xshock = 0.0;
        fctshk = 1.0;
    end

    if th~=0
        Re = Reth/th;
    else
        Re = Reth;
    end

    x = reshape(x, 1, []);
    Ue_tilde = reshape(Ue_tilde, 1, []);
    o = ones(size(x));
    A = [x; Ue_tilde; o; o];

    f = fopen(fpath, 'w');

    fprintf(f, '\t%6.4f\t%5.3e\t0.0000  \t\t\t |  Mach   Reyn   Rot\n', M, Re);
    fprintf(f, '\t0.00                            \t\t\t |  Hwrat\n');
    fprintf(f, '\t%6.4f\t20.0000                 \t\t\t |  Xtrip  Ncrit\n', xtrip);
    fprintf(f, '\t%8.6f\t%8.6f\t%8.6f    \t\t |  TH DS CT\n', th, ds, sqrt(ct));
    fprintf(f, '\t%8.6f\t%8.6f\t%8.6f    \t\t |  SCC SKP SKD\n', Kcorr, Kp, Kd);
    fprintf(f, '\t%8.6f\t%8.6f         \t\t\t |  XSHOCK FCTSHK\n', xshock, fctshk);
    fprintf(f, '\t%8.6f\t%8.6f\t%8.6f\t%8.6f\n', A);
    fclose(f)

end