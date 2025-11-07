function plot_entropy_budget(labels, varargin)

    groupLabels = {};
    for ir = 1:length(varargin)
        region = varargin{ir};
        stackData(ir, 1, 1) = region.e_conv;
        stackData(ir, 1, 2) = region.e_unst;
        stackData(ir, 1, 3:5) = 0;
    
        stackData(ir, 2, 1:2) = 0;
        stackData(ir, 2, 3) = region.e_diss;
        stackData(ir, 2, 4) = region.e_irrev;
        stackData(ir, 2, 5) = region.e_rev;
        groupLabels{ir} = labels(ir);
    end
    
    h = plotBarStackGroups(stackData, groupLabels);
    c = colororder;
    c = c(1:5,:);
    c = repelem(c,size(h,1),1); 
    c = mat2cell(c,ones(size(c,1),1),3);
    set(h,{'FaceColor'},c);
    legend(h(1,:),'\epsilon_S','\epsilon_{unst}','\epsilon_\phi','\epsilon_{irrev}','\epsilon_{rev}','Location','northeast')
    set(gca,'FontSize',12)
    

    set(gca, 'XTickLabelRotation',20)

end