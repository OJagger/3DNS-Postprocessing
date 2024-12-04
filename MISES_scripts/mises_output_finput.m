% mises output input variables (from a turbostream file)

function [f]=mises_output_finput(job,mode)

% mode = 1 constant DeHaller & load_coef
% mode = 2 constant streamtube contraction & load coef
% mode = 3 constant streamtube contraction & DH
% mode = 4 streamtube contraction set up to throat & DH

if exist('mode','var')==0
	mode = 1;
end

f.mode = mode;

% if job.Mrel_inlet < 1
% 	job.outname
% 	
% 	g=ts_read_hdf5([job.rjm_directory job.outname]);
% 	nj_mid = round(0.5*g{1}.attribute.nj);
% 	[bid_inlet,bid_outlet]=ts_inlet_bids(g);
% 	inlet=ts_structured_cut(g,bid_inlet,1,1,nj_mid-10,nj_mid+10,1,'en');
% 	outlet=ts_structured_cut(g,bid_outlet,'en','en',nj_mid-10,nj_mid+10,1,'en');
% 	
% 	inlet = ts_secondary(inlet);
% 	outlet = ts_secondary(outlet);
% 	
% 	if job.run_rotor == 1  % check what to do for stator
% 		f.AVDR = job.Ax;
% 	end
% 	
% 	
% 	[f.Pin,~]=ts_area_average(inlet,'P',3);
% 	[f.Pout,~]=ts_area_average(outlet,'P',3);
% 	
% 	if job.run_rotor == 1
% 		%     [f.M,~]=ts_mass_average(inlet,'M_rel',3);
% 		f.M = job.Mrel_inlet;
% 		%     [f.Alpha,~]=ts_mass_average(inlet,'Alpha_rel',3);
% 		f.Alpha = abs(job.alpha_rel_inlet);
% 		[f.V,~]=ts_mass_average(inlet,'V_rel',3);
% 		[f.Poin,~]=ts_mass_average(inlet,'Po_rel',3);
% 		[f.M2,~]=ts_mass_average(outlet,'M_rel',3);
% 		[f.Alpha2,~]=ts_mass_average(outlet,'Alpha_rel',3);
% 		[f.V2,~]=ts_mass_average(outlet,'V_rel',3);
% 	elseif job.run_rotor == 0
% 		%     [f.M,~]=ts_mass_average(inlet,'M',3);
% 		f.M = job.M_in;
% 		[f.Alpha,~]=ts_mass_average(inlet,'Alpha',3);
% 		f.Alpha = abs(job.alpha_inlet);
% 		[f.V,~]=ts_mass_average(inlet,'V',3);
% 		[f.Poin,~]=ts_mass_average(inlet,'Po',3);
% 		[f.M2,~]=ts_mass_average(outlet,'M',3);
% 		[f.Alpha2,~]=ts_mass_average(outlet,'Alpha',3);
% 		[f.V2,~]=ts_mass_average(outlet,'V',3);
% 	end
% elseif job.Mrel_inlet > 1
	
%% create intial guess based on initial dH
	ga=1.4;
	f.M = job.Mrel_inlet;
	f.Alpha = abs(job.alpha_rel_inlet);
	f.Alpha_abs = job.alpha_inlet;
	% just a reference not important
	Poin = 1./ga; % mises definition
	f.Poin = Poin;
	Pin = Poin.*(1+(ga-1)/2.*job.Mrel_inlet.^2).^(-ga/(ga-1));
	f.Pin = Pin;
	
	% assuming isentropic flow
	Poout  =  Poin;
	% assume polytropic efficiency of 0.95
	if isfield(job,'np')==0
		np = 0.95;
	else
		np = job.np;
	end
	f.np = np;
	% polytropic exponent
	n = np./(np-(ga-1)./ga);
	
	% 	if isfield(job,'DeHaller')==1;
	P_ratio = (1+(ga-1)./2.*f.M.^2.*(1-job.DeHaller.^2))^(ga.*np./(ga-1));
	job.P_ratio = P_ratio;
	f.P_ratio = P_ratio;
	Pout = job.P_ratio.*Pin;
	% 	elseif isfield(job,'P_ratio')==1;
	% 	Pout = job.P_ratio.*Pin;
	% 	end
	f.Pout = Pout;
	% calculate exit mach number
	syms Mout
	% 	Mout = vpasolve(Pout./Poout == (1+(ga-1)/2.*Mout.^2).^(-ga/(ga-1)),Mout,[0 1])
	V_cpTo1 = sqrt(ga-1).*f.M.*(1+(ga-1)./2.*f.M.^2).^(-0.5);
	
	%% note this assumes rel and absolute angles are of opposite sign
	M_theta = f.M*sind(f.Alpha) + f.M*cosd(f.Alpha)*tand(job.alpha_inlet);
	
	U_cpTo1 = sqrt(ga-1).*M_theta.*(1+(ga-1)./2.*M_theta.^2).^(-0.5);
	f.U_cpTo1 = U_cpTo1;
	
	flow_coefficient = V_cpTo1*cosd(f.Alpha)./U_cpTo1;
    
    %% calculate exit Mach number triangles based on mode
	
	if mode == 1
		Mout = vpasolve(job.DeHaller.*V_cpTo1 == sqrt(ga-1).*Mout.*(1+(ga-1)./2.*Mout.^2).^(-0.5),Mout,[0 1]);
		f.M2 = double(Mout);
		T_ratio = (1+(ga-1)/2.*f.M.^2)./(1+(ga-1)/2.*f.M2.^2);
		rho_ratio = job.P_ratio./T_ratio;
		rho_ratio = job.P_ratio.^(1./n);
		f.rho_ratio = rho_ratio;
		mass_nondim1 = ga/sqrt(ga-1)*f.M*(1+(ga-1)/2*f.M^2)^(-0.5*(ga+1)/(ga-1));
		mass_nondim2 = ga/sqrt(ga-1)*f.M2*(1+(ga-1)/2*f.M2^2)^(-0.5*(ga+1)/(ga-1));
		mass_nondim_P1 = ga/sqrt(ga-1)*f.M*(1+(ga-1)/2*f.M^2)^(0.5);
		mass_nondim_P2 = ga/sqrt(ga-1)*f.M2*(1+(ga-1)/2*f.M2^2)^(0.5);
		%% no contraction
		f.Alpha2 = acosd(mass_nondim_P1./mass_nondim_P2.*Pin./Pout.*cosd(f.Alpha));
		%% constant Mx
		f.Alpha2 = acosd(f.M.*cosd(f.Alpha)./f.M2);
		Ax_ratio = mass_nondim_P1./mass_nondim_P2.*Pin./Pout.*cosd(f.Alpha)./cosd(f.Alpha2);
		%% load coef, flow coefficient
		for n = 1:length(flow_coefficient)
			alpha_lim=-atand(job.load_coef./flow_coefficient(n)-tan(f.Alpha.*pi/180));
			syms alpha2
			alpha2=vpasolve( job.DeHaller.*cos(alpha2)./cos(f.Alpha.*pi/180) ...
				== (tan(f.Alpha.*pi/180)-job.load_coef/flow_coefficient(n))./tan(alpha2),alpha2,[(alpha_lim-5).*pi/180,pi/2-0.01]);
			
			Alpha2(n) = double(alpha2.*180./pi);
		end
		f.Alpha2 = Alpha2;
		Vx_ratio=job.DeHaller./cos(f.Alpha.*pi/180).*cosd(Alpha2);
		
		flow_coefficient.*(tand(f.Alpha)-Vx_ratio.*tand(Alpha2));
		alpha_target = abs(-atand(1/Vx_ratio*(-tand(f.Alpha)-job.load_coef/flow_coefficient)));
		% calculate Area ratio to achieve this
		Ax_ratio = 1./rho_ratio.*(1./Vx_ratio);
		f.load_coef = job.load_coef;
		f.throttle_phi = flow_coefficient;
		f.DeHaller = job.DeHaller;
		f.AVDR = Ax_ratio;
		f.Ax = Ax_ratio;
		
	elseif mode ==2
		
		
		DeHaller1 = job.DeHaller+0.1;
		DeHaller2 = job.DeHaller;
		while abs(DeHaller1 - DeHaller2)>0.001
			
			DeHaller1 = DeHaller2;
			p_ratio = (1+(ga-1)./2.*f.M.^2.*(1-DeHaller1.^2)).^(ga.*np./(ga-1));
			
			syms Mout
			Mout = vpasolve(DeHaller1.*V_cpTo1 == sqrt(ga-1).*Mout.*(1+(ga-1)./2.*Mout.^2).^(-0.5),Mout,[0 1]);
			f.M2 = double(Mout);
			
			rho_ratio = p_ratio.^(1./n);
			f.rho_ratio = rho_ratio;
			f.P_ratio = p_ratio;
			
			
			Ax_ratio = job.f.Ax;
			Vx_ratio = 1./Ax_ratio.*1./rho_ratio;
			alpha_lim=-atand(job.load_coef./flow_coefficient-tan(f.Alpha.*pi/180));
			syms alpha2
			alpha2=vpasolve( Vx_ratio ...
				== (tan(f.Alpha.*pi/180)-job.load_coef/flow_coefficient)./tan(alpha2),alpha2,[(alpha_lim-5).*pi/180,pi/2-0.01]);
			Alpha2 = double(alpha2.*180./pi);
			
			f.Alpha2 = Alpha2;
			job.DeHaller=Vx_ratio.*cos(f.Alpha.*pi/180)./cosd(Alpha2);
			DeHaller2 = job.DeHaller;
			alpha_target = abs(-atand(1/Vx_ratio*(-tand(f.Alpha)-job.load_coef/flow_coefficient)));
			f.load_coef = job.load_coef;
			f.throttle_phi = flow_coefficient;
			f.DeHaller = job.DeHaller;
			f.AVDR = Ax_ratio;
			f.Ax = Ax_ratio;
		end
		
	elseif mode == 3
		
		Mout = vpasolve(job.DeHaller.*V_cpTo1 == sqrt(ga-1).*Mout.*(1+(ga-1)./2.*Mout.^2).^(-0.5),Mout,[0 1]);
		f.M2 = double(Mout);
		T_ratio = (1+(ga-1)/2.*f.M.^2)./(1+(ga-1)/2.*f.M2.^2);
		rho_ratio = job.P_ratio./T_ratio;
		rho_ratio = job.P_ratio.^(1./n);
		f.rho_ratio = rho_ratio;
		
		Ax_ratio = job.f.Ax;
		Vx_ratio = 1./Ax_ratio.*1./rho_ratio;
		alpha2 = acos(Vx_ratio.*cos(f.Alpha.*pi/180)./job.DeHaller);
		Alpha2 = double(alpha2.*180./pi);
		f.Alpha2 = Alpha2;
			
		%% load coef, flow coefficient
		for n = 1:length(flow_coefficient)
			
			syms load_coef
			load_coef=vpasolve( job.DeHaller.*cos(alpha2)./cos(f.Alpha.*pi/180) ...
				== (tan(f.Alpha.*pi/180)-load_coef/flow_coefficient(n))./tan(alpha2),load_coef,[0.25,0.5]);
			
		end
		
		job.load_coef = double(load_coef);
		
		% calculate Area ratio to achieve this
		f.load_coef = job.load_coef;
		f.throttle_phi = flow_coefficient;
		f.DeHaller = job.DeHaller;
		f.AVDR = Ax_ratio;
		f.Ax = Ax_ratio;
        
    elseif mode == 4
        
        Mout = vpasolve(job.DeHaller.*V_cpTo1 == sqrt(ga-1).*Mout.*(1+(ga-1)./2.*Mout.^2).^(-0.5),Mout,[0 1]);
        f.M2 = double(Mout);
        T_ratio = (1+(ga-1)/2.*f.M.^2)./(1+(ga-1)/2.*f.M2.^2);
        rho_ratio = job.P_ratio./T_ratio;
        rho_ratio = job.P_ratio.^(1./n);
        f.rho_ratio = rho_ratio;
        
        %% set to one!
        Ax_ratio = 1;
        %%
        Vx_ratio = 1./Ax_ratio.*1./rho_ratio;
        alpha2 = acos(Vx_ratio.*cos(f.Alpha.*pi/180)./job.DeHaller);
        Alpha2 = double(alpha2.*180./pi);
        f.Alpha2 = Alpha2;
        
        %% load coef, flow coefficient
        for n = 1:length(flow_coefficient)
            
            syms load_coef
            load_coef=vpasolve( job.DeHaller.*cos(alpha2)./cos(f.Alpha.*pi/180) ...
                == (tan(f.Alpha.*pi/180)-load_coef/flow_coefficient(n))./tan(alpha2),load_coef,[0.25,0.5]);
            
        end
        
        job.load_coef = double(load_coef);
        
        % calculate Area ratio to achieve this
        f.load_coef = job.load_coef;
        f.throttle_phi = flow_coefficient;
        f.DeHaller = job.DeHaller;
        f.AVDR = job.f.Ax;
        f.Ax = job.f.Ax;
        
        
    end
% end


% calculate Vta01
ga=1.4;
f.Vtao1 = f.M/(sqrt(1+(ga-1)/2*f.M^2))*tand(abs(f.Alpha))/(sqrt(1+tand(abs(f.Alpha))^2));

% calcualte Reynolds number
% mou = g{1}.av.viscosity;
% [f.ro,~]=ts_mass_average(inlet,'ro',3)
% f.Re = f.ro*f.V*job.chord/mou;
f.Re = 1*10^6;

end