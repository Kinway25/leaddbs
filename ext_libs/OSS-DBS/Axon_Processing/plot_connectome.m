tracts = dir_without_dots('/media/interscan/BackupKB/JS_VTAs_SVDs/sub-8qsELGyn148nzUCfX2r4Tk/SVD/PetersenUD4_WMH');

for tract_i = 1:size(tracts,1)

    ea_plot_streamtubes([tracts(tract_i).folder,filesep,tracts(tract_i).name]);
    % purple
    %set(mytract,'FaceColor',[0.4940 0.1840 0.5560],'FaceAlpha',1.0,'EdgeColor','none')
end
