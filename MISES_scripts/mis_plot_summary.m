% function to plot overall summary plots after mis_store_AtA1 has been run

function mis_plot_summary(directory,var_x,var_y,fold)
% inputs: datum file you want to look at
% var_y: mis store variables to plot in y_axis
% var_x: mis_store variable to plot in x_axis
% fold: pick particular folders

%% linestyle, markerstyle and colors
color = [0,0,1;1,0,0;0,0.5,0;1,0,1;0.5,0,0.5;0.5,0.5,0.5;0,0,0;0.25,0.25,0.25;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0];
Linestyle = {'-','--',':','-.','-','-','-','-','-','-'};
markerstyle = {'x','o','+','*','s','d','^','v','>','<','x','o','+','*','s','d','^','v','>','<'};

%% define variables to plot on y-axis and x-axis if not pre-defined in inputs to function

if exist('var_y','var')==0 || isempty(var_y) == 1
    var_y = {'dPo_Po','M_max','P','P_is'};
end

if exist('var_x','var')==0 || isempty(var_x) == 1
    var_x = {'AtA1'};
end

if exist('fold','var') == 0 
   fold = []; 
end


%% find directory files
dir_dat =  '/mnt/Disk2/dl467/Documents/MISES_summary/';
A = [dir([dir_dat directory(1:3) '/' directory '/'])];
F = cell(length(A),1); for n = 1:length(A); F{n} = A(n).name; end;
% remove rogue directories based on dots
F(cellfun(@isempty,strfind(F,['.'])) == 0) = [];
if isempty(fold) == 0
    for kk = 1:length(fold)
        F_temp{kk} = F(cellfun(@isempty,strfind(F,[fold{kk}])) == 0);
        if length(F_temp{kk})~=1
            F_temp{kk} = F_temp{kk}(1); % this is to distinguish between 'Ax' and 'Axt'
        end
    end
    F = F_temp;
end

%% create number of figures based on number of variables to be plotted
for k = 1:length(var_y)
    h(k) = figure(); grid on; hold on;
    xlabel(var_x{1}); ylabel(var_y{k}); title(directory);
end

%% plot mis store variables
nn = 1;
out_dat=str2double(regexp(directory,'[\d.]+','match'));

for i = 1:length(F)
    B = dir([dir_dat directory(1:3) '/' directory '/' char(F{i})]);
    F2 = cell(length(B),1);
    for n = 1:length(B); F2{n} = B(n).name; end;
    F2(cellfun(@isempty,strfind(F2,['.mat'])) == 1) = [];
    for j = 1:length(F2)
        out=str2double(regexp(F2{j},'[\d.]+','match'));
        for nj=1:length(out_dat)
            if out(nj) ~= out_dat(nj)
                nk = nj;
                break
            end
        end
        
        if length(out) == 6
            def = {'M' ; 'A' ; 'SC' ; 'TC' ; 'PSI' ; 'A'};
        else
            def = {'M' ; 'A' ; 'SC' ; 'TC' ; 'PSI' ; 'Ax' ; 'A'};
        end
        leg{nn} = [def{nk} num2str(out(nk))];
        nn = nn+1;
         load([dir_dat directory(1:3) '/' directory '/' char(F{i}) '/' char(F2{j})])
        for k = 1:length(var_y)
            for l = 1:length(var_x)
            figure(h(k));
            plot(mis.(var_x{l}),mis.(var_y{k}),markerstyle{nn})
            legend(leg)
            end
        end
    end
end

end
        


