function ea_plot_streamtubes(varargin)

    % example 
    % ea_plot_streamtubes('/home/forel/Documents/GitHub/leaddbs/connectomes/dMRI_MultiTract/AysuFilteredByROIDownsampledBy1/GPe_STN_caudal_250_smooth_right.mat')

    fiberActivationProb = varargin{1};
    load(fiberActivationProb);

    if nargin >=2
        numfibers = varargin{2};
    else
        numfibers = size(idx,1);
    end

    if numfibers > 500
        downsamplefactor = 5;
    else
        downsamplefactor = 1;
    end

    col = [1,0,0];

    %% downsampling
    fibersnew=mat2cell(fibers(:,1:3),idx);
    fibersnew = cellfun(@(f,len) f(round(linspace(1,len,round(len/downsamplefactor))),:), fibersnew, num2cell(cellfun(@(p) size(p,1), fibersnew)), 'UniformOutput', 0);

    for fiber_i = 1:size(fibersnew,1)
        mytract = streamtube(fibersnew(fiber_i),0.1);
        %set(mytract,'FaceColor',[0.83 0.25 0.25],'FaceAlpha',1.0,'EdgeColor','none')
        
        %set(mytract,'FaceColor',[0.4940 0.1840 0.5560],'FaceAlpha',1.0,'EdgeColor','none')
        set(mytract,'FaceColor',[0.3940 0.0840 0.4560],'FaceAlpha',1.0,'EdgeColor','none')
        %set(mytract,'FaceColor',[0.6940 0.3840 0.7560],'FaceAlpha',1.0,'EdgeColor','none')
    end

    
