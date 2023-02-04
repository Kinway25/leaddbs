function spatio_corr_mat = fiber_spatial_rev3(fibers,trgt_coor,sig,vcnty_thr)
%% Fiber Spatial Analysis: Rev 3 
% Min Jae Kim (mkim@bwh.harvard.edu)
% Last Edit: 01/30/2023
%% Loading Fibers
%load('gpe2stn_sm_right.mat');

%% Organizing Spatial Coordinates of Each Fibers 
num_fibers=length(unique(fibers(:,4))); % number of fibers
fiber_idxi=cell(num_fibers,2); % data structure
for i=1:num_fibers
    fiber_idxi{i,1}=i; %1st column: fiber number 
    idxi=find(fibers(:,4)==i);
    coors=fibers(idxi,:);
    fiber_idxi{i,2}=coors(:,1:3); %2nd column: fiber coordinates
end

%% Step 0: Initializing Parameters 
method='Spearman'; % Method for Correlation: (1) Spearman or (2) Pearson
MNI_temp=ea_load_nii([ea_space,'t1.nii']);
dim=MNI_temp.dim;
%% Method 1: Spatial Correlation Method
spatio_corr_mat=zeros(num_fibers,num_fibers);
%spatio_corr_mat=zeros(20,20);% sample matrix for first 20 fibers
%% Processing Fibers (3D Gaussianization)
fib_gaus_cell=cell(num_fibers,2);
for i=1:num_fibers
    smp1=fiber_idxi{i,2};
    [smp1_gauss,idxi_1]=fib_gauss(smp1,sig,trgt_coor,vcnty_thr,dim);
    fib_gaus_cell{i,1}=smp1_gauss;
    fib_gaus_cell{i,2}=idxi_1;
    disp("Processing Fiber Number #"+string(i))
end 

%% Performing Spatial Correlation (3D Gaussianization)
for i=1:num_fibers
    smp1=fib_gaus_cell{i,1};
    idxi_1=fib_gaus_cell{i,2};
    tic
    for j=1:num_fibers
        smp2=fib_gaus_cell{j,1};
        idxi_2=fib_gaus_cell{j,2};
        [C,ia,ib] = union(idxi_1,idxi_2);
        [C_1,i_1,~]=intersect(C,idxi_1);
        [C_2,i_2,~]=intersect(C,idxi_2);
        
        temp1=zeros(length(C),1);
        temp1(i_1)=smp1;
        temp2=zeros(length(C),1);
        temp2(i_2)=smp2;
        
        rho=corr(temp1,temp2,'Type',method,'Rows','complete');
        spatio_corr_mat(i,j)=rho;
    end
    toc
end
%ea_spatial_corr
% Visualization
figure
heatmap(spatio_corr_mat);
grid off
colormap jet
xlabel("Fiber Number")
ylabel("Fiber Number");
title("Pairwise Spatial Correlation (r) Between Fibers");
end
% %% Method 2: Eucledian Distance Method
% 
% tic
% euc_dist_mat=zeros(num_fibers,num_fibers);
% for i=1:num_fibers
%     smp1=fiber_idxi{i,2};
%     for j=1:num_fibers
% 
%         if i == 41 && j == 40
%             disp("check")
%         end
% 
%         smp2=fiber_idxi{j,2};
%         mean_euc=p_euc(smp1,smp2,trgt_coor,vcnty_thr);
%         euc_dist_mat(i,j)=mean_euc;
% 
%         % take the larger distance 
%         if j < i
%             if euc_dist_mat(i,j) < euc_dist_mat(j,i)
%                 euc_dist_mat(i,j) = euc_dist_mat(j,i);
%             else
%                 euc_dist_mat(j,i) = euc_dist_mat(i,j);
%             end
%         end
% 
%     end 
% end 
% toc
% 
% % Visualization
% figure
% heatmap(euc_dist_mat);
% grid off
% colormap jet
% xlabel("Fiber Number")
% ylabel("Fiber Number");
% title("Mean Pairwise Eucledian Distance (mm) Between Fibers");
% %% 