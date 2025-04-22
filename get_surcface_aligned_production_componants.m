function [Pr_s, Pr_n, Pr_t] = get_surcface_aligned_production_componants(s)

    % DUDX = s.oGridProp('DUDX');
    % DUDY = s.oGridProp('DUDY');
    % DVDX = s.oGridProp('DVDX');
    % DVDY = s.oGridProp('DVDY');
    % 
    % roUddUdd = cas_trip.meanFlow.oGridProp('roUddUdd');
    % roVddVdd = cas_trip.meanFlow.oGridProp('roVddVdd');
    % roUddVdd = cas_trip.meanFlow.oGridProp('roUddVdd');
    % 
    % pr_tensor = zeros(size(roUddUdd));
    % pr_tensor = repmat(pr_tensor,[1 1 3 3]);
    % S_tensor = zeros(size(pr_tensor));
    % tau_tensor = zeros(size(pr_tensor));
    % 
    % tau_tensor(:,:,1,1) = -roUddUdd;
    % tau_tensor(:,:,2,2) = -roVddVdd;
    % tau_tensor(:,:,1,2) = -roUddVdd;
    % tau_tensor(:,:,2,1) = tau_tensor(:,:,1,2);
    % 
    % S_tensor(:,:,1,1) = DUDX;
    % S_tensor(:,:,1,2) = 0.5*(DUDY + DVDX);
    % S_tensor(:,:,2,1) = S_tensor(:,:,1,2);
    % S_tensor(:,:,2,2) = DVDY;

    for i=1:size(DUDX,1)
        n = s.n(:,i);


    end



end