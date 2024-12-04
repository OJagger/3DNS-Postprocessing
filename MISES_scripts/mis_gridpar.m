function mis_gridpar(directory,n_type)
% create gridpar file quickly when re-running
%n_type= 1 : n=100;
%n_type= 2 : n=75;
%n_type= 3 : D_scl = 0.07;
%n_type= 4 : n=125;

% Write gridpar file
fid = fopen([directory 'gridpar.mises'],'w');
fprintf(fid,'%s\n','T T');
fprintf(fid,'%s\n','60   45');
fprintf(fid,'%s\n','24');
fprintf(fid,'%s\n','1.500000');
fprintf(fid,'%s\n','0.800000');
if n_type == 1
	n=100;
	fprintf(fid,'%s\n',[num2str(n-1) '    0.10000    0.900000    1.000000']);
elseif n_type == 2
	n=75;
	fprintf(fid,'%s\n',[num2str(n-1) '    0.10000    0.900000    1.000000']);
elseif n_type == 3
	n=100;
	fprintf(fid,'%s\n',[num2str(n-1) '    0.07000    0.900000    1.000000']);
elseif n_type == 4
    n=125;
	fprintf(fid,'%s\n',[num2str(n-1) '    0.10000    0.900000    1.000000']);
elseif n_type == 5
    n=75;
	fprintf(fid,'%s\n',[num2str(n-1) '    0.0700    0.900000    1.000000']);	
elseif n_type == 6
    n=125;
	fprintf(fid,'%s\n',[num2str(n-1) '    0.0500    0.900000    1.000000']);	
elseif n_type == 7
    n=85;
	fprintf(fid,'%s\n',[num2str(n-1) '    0.0850    0.900000    1.000000']);
elseif n_type == 8
	n=100;
	fprintf(fid,'%s\n',[num2str(n-1) '    0.08500    0.900000    1.000000']);	
end
fprintf(fid,'%s\n','1.000000    1.000000    0.000000    1.000000    1.000000    0.000000');
fclose(fid);

end