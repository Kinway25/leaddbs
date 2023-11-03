function mytract = till_addtract(myfile,thickness,col,alph,downsamplefactor,numfibers)
    load(myfile);
    %% old and slow
%     fibs = unique(fibers(:,4));
%     for k = 1:length(fibs)
%        fibersnew{k,1} =  fibers(fibers(:,4) == fibs(k),1:3);
%     end    
    %% new and fast
    fibersnew=mat2cell(fibers(:,1:3),idx);
    %% downsampling
    fibersnew = cellfun(@(f,len) f(round(linspace(1,len,round(len/downsamplefactor))),:), fibersnew, num2cell(cellfun(@(p) size(p,1), fibersnew)), 'UniformOutput', 0);
    %% reduce number of visualized fibers
    if ~isempty(numfibers) && length(fibersnew) > numfibers
        myfibs = ceil(linspace(1, length(fibersnew), numfibers));        
    else
        myfibs = 1:length(fibersnew);
    end
    %%    
    mytract = streamtube(fibersnew(myfibs),thickness);
    set(mytract,'FaceColor',col,'FaceAlpha',alph,'EdgeColor','none')
end