function [n_type]=mis_run_inner(directory,kspec,n_type)
% function to run inner iteration

cd(directory)

if exist('n_type','var')==0 || isempty(n_type)==1
	n_type = 1;
end

n_type_orig = n_type
mis_gridpar(directory,n_type_orig);

n_type_lim = 7;

if exist('kspec','var')==0 || isempty(kspec)==1
	kspec = 1;
end

if n_type_orig == 1
if kspec == 1
	[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle.py"','-echo');
elseif kspec == 17
	[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle2.py"','-echo');
end
else % if you have to run again (say second polar in mis_const_speed)
	[~,~]=system(['rm ' directory 'polarx.mises']);
	[~,~]=system(['rm ' directory 'idat.mises']);
end

% Read in flow file and re-run if not converged
if exist([directory 'polarx.mises'],'file') ~= 0
	[Polarx, Ises] = mis_read_polarx('mises',directory);
	if isstruct(Polarx) ~= 1
		disp('*******Run again*********');
	    [~,~]=system(['rm ' directory 'polarx.mises']);
		while isstruct(Polarx) ~= 1 && n_type < n_type_lim
			n_type = n_type + 1;
			disp(['****Create new gridpar ' num2str(n_type) '****']);
			mis_gridpar(directory,n_type);
			if kspec == 1
				[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle.py"','-echo');
			elseif kspec == 17
				[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle2.py"','-echo');
			end
			[Polarx, Ises] = mis_read_polarx('mises',directory);
		end
		if n_type>n_type_lim && isstruct(Polarx) ~= 1
			disp('*******Have to run manually*********');
			while isstruct(Polarx) ~= 1
				[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle_manual.py"','-echo');
				[Polarx, Ises] = mis_read_polarx('mises',directory);
			end
		end
% 		n_type = n_type_orig; mis_gridpar(directory,n_type);
	else
		n_type = 2;
	end
else
	disp('File Not Found')
	disp('*******Run again*********');
	Polarx = [];
	while isstruct(Polarx) ~= 1 && n_type < n_type_lim
		n_type = n_type + 1;
		disp(['****Create new gridpar ' num2str(n_type) '****']);
		mis_gridpar(directory,n_type);
		if kspec == 1
			[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle.py"','-echo');
		elseif kspec == 17
			[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle2.py"','-echo');
		end
		[Polarx, Ises] = mis_read_polarx('mises',directory);
	end
	if n_type==n_type_lim && isstruct(Polarx) ~= 1
		disp('*******Have to run manually*********');
		while isstruct(Polarx) ~= 1
			if kspec == 1
				[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle_manual.py"','-echo');
			elseif kspec == 17
				[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle_manual2.py"','-echo');
			end
			[Polarx, Ises] = mis_read_polarx('mises',directory);
		end
	end
% 	n_type = n_type_orig; mis_gridpar(directory,n_type);
	p = [];
	% 			return
end

end