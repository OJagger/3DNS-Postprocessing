function cases = meanflow_convergence(basecase, runs, prop)

    if ~iscell(basecase)
        cases = {};
    

        for run=runs
            casenow = feval( class(basecase), basecase.casepath, run);
            casenow.readMeanFlow;
            cases{end+1} = casenow;
        end
    else
        cases = basecase;
    end

    q = [];
    for run = 1:length(cases)
        q(run) = cases{run}.meanFlow.(prop);
    end
    
    figure
    plot(runs, q, '-ks', 'MarkerFaceColor','k');
    grid on
    xlabel('Run')
    ylabel(prop);


end