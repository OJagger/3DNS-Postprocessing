function slice_animation(casename,run,var,nworkers)
var = string(var);
fprintf('Processing case %s, run %d, %s\n',casename,run, var)
p = gcp('nocreate');
if isempty(p)
    parpool(nworkers);
end
hfcase = DNS_case(casename, run);
%hfcase = DNS_case('cwl90_window_turb_clean',3);
%runs =[5];
%hfcase.readKSlices([],61:85);
imgfolder = fullfile(hfcase.casepath,sprintf('run%d',run),'animation_images',var);
%%
if ~exist(imgfolder, 'dir')
       mkdir(imgfolder);
end

%%
parfor i=1:hfcase.nSlices
    fprintf('Plotting slice %d/%d\n',i,hfcase.nSlices)
    if ~exist(fullfile(imgfolder,sprintf('img_%03d.png',i)),'file')
    if strcmp(var,"k_Mach")
        slice = hfcase.readSingleKSlice(i);
        slice2kPlot(slice, hfcase.blk, 'M', fullfile(imgfolder,sprintf('img_%03d.png',i)), [0 1.4], 'M');
    elseif strcmp(var,"j_tau")
		slice = hfcase.readSingleJSlice(i);
		slice2jPlot(slice, 'tau_w', fullfile(imgfolder, sprintf('img_%03d.png',i)), [0 800], '\tau_w');
    elseif strcmp(var,"bl_vort")
        slice = hfcase.readSingleKSlice(i)
        slice2kPlot_closeup(slice, hfcase.blk, 'vortZ', fullfile(imgfolder, sprintf('img_%03d.png',i)), 1e5*[-0.7 0.7], '\omega_z')
        
    end
    end
    %clear slice;
end

end

