function v = ismax(prof)

    v(1:length(prof)) = false;
    d = diff(prof);
    s = d(2:end) .* d(1:end-1);
    v(2:end-1) = s<0 & d(1:end-1) > 0;

end