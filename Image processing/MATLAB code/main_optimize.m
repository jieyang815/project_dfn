clc
clear
close all

% -------------------------------------------------------------------------
% A bi-directional evaluation is performed on the fiber segmentation results
% three parameters are evaluated: number of clusters (kmeans, for binarization)
% number of tiles and clip limit (CLAHE, for enhancement)
% -------------------------------------------------------------------------

%% --------  Preliminary setup -------------------------------------
resize_scale = 100; % scaling up factor used in manual segmentation

% ---- Binarization ---------------------------------------------------------
seg_type = 2; % parameter to switch binarization type
% 1: Otsu's method
% 2: K-means clustering

% ------ Denoise ----------------------------------------------------------
size_close = 2; % close small holes after binarization

% ---- Postprocessing: connect broken fibers -------------------------------------------
error_acc = 15; % slope tolerance (°): criteria to select fibers with similar slope
error_ll = 15; % parallel slope tolerance (°) : criteria to define if fibers are parallel

% ----- Skeletonization and nodes detect -------------------------------
min_quality = 0.25; % for Harris corner detection: minimum accepted quality of corner metric value

% --- Evaluation ----------------------------------------------------------
epsilon = 0.01; % sampling interval
limit_Nmax = 0.25; % limit for negative metric, to suppress false negative

%% -- Import data: original image and ground truth -----------------------
[data_img,~] = readPath('Select Sample data folder',resize_scale,epsilon);

%% ---- Select folders to save reconstruction results --------------------
path_eva = uigetdir(' ','Select the folder to save evaluation results');
path_DFN = uigetdir(' ','Select the folder to save the geometry of reconstructed fiber network');

%% Parameter analysis: CLAHE and kmeans
% -------------------------------------------------------------------------
% CLAHE is used to enhance the image by increasing the contrast between objects
% and background. Clip limit and number of tiles are to be determined
% -------------------------------------------------------------------------
% Kmeans clustering is used to binarize the image, and the number of clusters
% need to be determined
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
tic
% -------------------------------------------------------------------------
% ---- Parameter to evaluate ----------------------------------------------
list_numcluster = 2:7; % number of clusters in kmeans clustering
list_numtile = [4 8 16 32 64];
list_cliplimit = 0.01:(1-0.01)/(length(list_numtile) - 1):1;

% ---- Visualize evaluation results -----------------------------------
f1 = figure('Name','FNR','NumberTitle','off');
f2 = figure('Name','FPR','NumberTitle','off');
% ---------------------------------------------------------------------------

% ---- Parameter analysis on kmeans and CLAHE begins ----------------------
for m = 1:length(list_numcluster)
    num_cluster = list_numcluster(m); % number of clusters to evaluate for kmeans
    
    fprintf(1, 'Evaluation for %.f clusters starts\n', num_cluster);
    
    % --- Initialize evaluation metric array --------------------------------
    list_Pmean = zeros(size(list_cliplimit,2),size(list_numtile,2),length(data_img)); % list of postive metric
    list_Nmean = list_Pmean; % list of negative metric
    dim = ndims(list_Pmean);
    for i = 1:length(data_img) % iterate through all images to find the average optimal value
        
        % --- import ground truth of image i ---------------------------
        xy_GT = data_img(i).GT;
        
        % --- import original image i ----------------------------------
        img_original = data_img(i).img;
        name_tosave = data_img(i).imgname;
        img_h = data_img(i).height;
        len_short = data_img(i).minlen;
        sigma = data_img(i).sigma;
        spacing = data_img(i).spacing;
        
        % --- Iteration of number of clusters to evaluate segmentation performance ------
        for j = 1:size(list_cliplimit,2)
            limit_clip = list_cliplimit(j);
            
            for k = 1:size(list_numtile,2)
                num_tile = list_numtile(k);
                
                % --- Segment the network from image ----------------------
                parameter_pre = [limit_clip num_tile]; % parameter for image enhancement
                xy_SEG = ...
                    FiberSeg(img_original,parameter_pre,num_cluster,...
                    size_close,len_short,min_quality,error_acc,error_ll);
                
                % ----- Perform bi-directional evaluation on segmentation results ----
                [Pmean,Nmean,~,~,~,~] = getMetric(xy_GT,xy_SEG,sigma,spacing);
                list_Pmean(j,k,i) = Pmean;
                list_Nmean(j,k,i) = Nmean;
                fprintf(1, 'Analysis for %s is completed with clip = %.2f, tiles = %.f\n', name_tosave,limit_clip,num_tile);
                
            end
            
        end
        
    end
    % --- Data processing -----------------------------------------------------
    [P_mean,N_mean,Pmax_temp,Pstd_temp,clipmax_temp,tilemax_temp] = ...
        dataProcess(list_Pmean,dim,list_cliplimit,list_numtile,list_Nmean,limit_Nmax);
    
    % ---- save parameter optimization result: average across all images -------------------------------------------
    eva_fullpath = fullfile(path_eva,['Eva_kNum',num2str(num_cluster,'%.f')]);
    save(eva_fullpath,'P_mean','N_mean')
    % --- plot: average positive metric of 5 images -----------------------------------------------
    FNR = 1 - P_mean;
    FPR = N_mean;
    
    figure(f1)
    subplot(ceil(length(list_numcluster)/2),2,m)
    evaluatePlot(FNR,list_numtile,list_cliplimit,num_cluster,size(data_img,2),0.1,0.4)
    
    % --- average negative metric of 5 images -----------------------------------------
    figure(f2)
    subplot(ceil(length(list_numcluster)/2),2,m)
    evaluatePlot(FPR,list_numtile,list_cliplimit,num_cluster,size(data_img,2),0.1,0.4)
    % --------------------------------------------------------------------------------------
    
    if m == 1 % first round of cluster evaluation
        Pmean_optimal = Pmax_temp;
        clip_optimal = clipmax_temp;
        numtile_optimal = tilemax_temp;
        cluster_optimal = num_cluster;
    else
        if Pmax_temp > Pmean_optimal % higher positive metric
            Pmean_optimal = Pmax_temp;
            clip_optimal = clipmax_temp;
            numtile_optimal = tilemax_temp;
            cluster_optimal = num_cluster;
        end
    end
    
    fprintf(1, 'Evaluation for %.f clusters finished\n', num_cluster);
    
end
% -------------------------------------------------------------------------
toc
% -------------------------------------------------------------------------

%% ---- Plot optimal results ----------------------------
% clip_optimal = 0.05;
% numtile_optimal = 32;
% cluster_optimal = 3;
type_pre = 2;
for i = 1:length(data_img)
    % --- import ground truth of image i ---------------------------
    xy_GT = data_img(i).GT;
    
    % --- import original image i ----------------------------------
    img_original = data_img(i).img;
    name_tosave = data_img(i).imgname;
    img_h = data_img(i).height;
    len_short = data_img(i).minlen;
    sigma = data_img(i).sigma;
    spacing = data_img(i).spacing;
    
    % --- Segment the network from image ----------------------
    parameter_optimal = [clip_optimal numtile_optimal];
    xy_SEG = ...
        FiberSeg(img_original,parameter_optimal,cluster_optimal,...
        size_close,len_short,min_quality,error_acc,error_ll);
    
    % ----- Perform bi-directional evaluation on segmentation results ----
    [Pmean_pmax,Nmean_pmax,Pmetric_pmax,Nmetric_pmax,gridGT_pmax,gridSEG_pmax]...
        = getMetric(xy_GT,xy_SEG,sigma,spacing);
    
    % ---- export reconstruction results ------------------------------------
    save([path_DFN,'\Data_reNetwork_img',num2str(i,'%.f')],'xy_SEG','clip_optimal','numtile_optimal','cluster_optimal')
    
    % ----- Visualize evaluation metric locally by mapping it to the network --------------------------------------
    f1 = figure('Name',name_tosave,'NumberTitle','off');
    figure(f1)
    t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
    ax = gobjects(1,3);
    
    % --- show original image ----------
    ax(1) = nexttile();
    imshow(img_original,[])
    title(name_tosave)
    axis equal
    
    % ---- map positive metric to the ground truth network ---------------
    ax(2) = nexttile();
    title_Pmean = ['FNR = ',num2str(1 - Pmean_pmax,'%.2f')];
    plotLocal(xy_SEG,gridGT_pmax,Pmetric_pmax,sigma,img_h,title_Pmean)
    axis equal
    
    % ---- map negative metric to the segmented network ------------------
    ax(3) = nexttile();
    title_Nmean = ['FPR = ',num2str(Nmean_pmax,'%.2f')];
    plotLocal(xy_GT,gridSEG_pmax,Nmetric_pmax,sigma,img_h,title_Nmean)
    axis equal
    
    % ---- set colormap --------------------------------------------------
    colormap(ax(1),'gray')
    
    fprintf(1, 'Plot for %s is completed\n', name_tosave);
    
    % --- save results as .fig & .png ------------------------------------
    savepath_full = fullfile(path_eva,['Evaluation_K_CLAHE_',name_tosave]);
    set(f1, 'Position', get(0, 'Screensize'));
    frame1    = getframe(f1);
    imwrite(frame1.cdata, [savepath_full,'.png'], 'png')
    %      print(f1,savepath_full,'-dpng')
    savefig([savepath_full,'.fig'])
end

% -------------------------------------------------------------------------
toc
% -------------------------------------------------------------------------

%% -------------- Function ---------------------------------------------------------
function [data_img,path_tosave] = readPath(text,resize_scale,epsilon)
% Input:
% text: command to select specified folder
% resize_scale: scaling factor during manual segmentation
% epsilon: sampling interval for interpolation

% Output:
% data_img: store the information of each image: ground truth, evaluation parameter etc.
% path_tosave: directory to save the results
path_tosave = uigetdir(' ',text); % select the data folder
% ----- original image path -------------------------
pattern_temp = fullfile(path_tosave,'*.tif');
path_original = dir(pattern_temp);

% ------ Ground truth path list---------------------------
pattern_temp = fullfile(path_tosave,'*.txt');
path_GT = dir(pattern_temp);

% ---- Read data -------------------------------
for i = 1:length(path_GT)
    % ---- import ground truth --------------------------------------------
    name_temp = path_GT(i).name;
    path_temp = fullfile(path_GT(i).folder, name_temp); % construct complete file path
    xy_GT = importTXT(path_temp,resize_scale);
    
    % ---- import original image ------------------------------------------
    name_temp = path_original(i).name;
    path_temp = fullfile(path_original(i).folder, name_temp); % construct complete file path
    [img_original,name_tosave,img_h,len_short,sigma,spacing] = ...
        importImg(name_temp,path_temp,epsilon);
    
    % ---- save image data ------------------------------------------------
    data_img(i).GT = xy_GT; % record coordinate of fibers in GT
    data_img(i).img = img_original;
    data_img(i).imgname = name_tosave;
    data_img(i).height = img_h; % record coordinate of fibers in GT
    data_img(i).minlen = len_short; % record coordinate of fibers in GT
    data_img(i).sigma = sigma; % record coordinate of fibers in GT
    data_img(i).spacing = spacing; % record coordinate of fibers in GT
end

end

function xy_GT = importTXT(path_temp,resize_scale)
file_txt = fopen(path_temp,'rt');
xy_GT = {};
while ~feof(file_txt) % iterate along each line in the file
    row = fgetl(file_txt);
    x_vector = str2num(row) ./ resize_scale; % resize x vector
    row = fgetl(file_txt);
    y_vector = str2num(row)  ./ resize_scale;
    xy_GT = [xy_GT;[y_vector' x_vector']];
end
end

function [img_original,name_tosave,img_h,len_short,sigma,spacing] = ...
    importImg(name_temp,path_temp,epsilon)

% ------ import original image --------------------------------------------
img_original = imread(path_temp);

% ---- construct file name to save data -------------------------------
name_cut = split(char(name_temp),'_'); % cut before file number for file name reconstruction
name_tosave = char(join(name_cut(1)));

% --- parameter depending on image size -------------------------------
[img_w,img_h] = size(img_original);
len_short = ceil(max(img_w,img_h) * 0.05); % minimum branch length on the skeleton, branch shorter than this limit will be removed
sigma = 0.025 * (max(img_w,img_h)); % calculate sensitive parameter for evaluation: set as 0.025*image size
spacing = sigma * epsilon; % sampling internal
end

function [P_mean,N_mean,max_Pmean,std_Pmean,var1_optimal,var2_optimal] = dataProcess(list_Pmean,dim,var1,var2,list_Nmean,limit_Nmax)
% Input;
% list_Pmean, list_Nmean: evaluation metrics of the images
% dim: dimension to conduct averaging
% var1, var2: parameters to evaluate
% limit_Nmax: maximum tolerance for false positive rate

% Output:
% P_mean, N_mean: average across all images
% max-Pmean: maximum positive metric
% std_Pmean: standard deviation of positive metric
% var1_optimal, var2_optimal: optimal parameter set that yield the best result

% ---- filter out negative metric > limit ------------------------------------
N_mean = mean(list_Nmean, dim); % average data from all the images
index_Nmax = N_mean < limit_Nmax;

% --- find maximum positive metric after filering ----------------------------------------
P_mean = mean(list_Pmean, dim); % average data from all the images
Pmean_filtered = P_mean .* index_Nmax;
std_Pmean = std(list_Pmean,0,dim); % standard deviation of dataset
max_Pmean = max(max(Pmean_filtered)); % maximum metric value
[row_max,col_max] = find(Pmean_filtered == max_Pmean,1);
var1_optimal = var1(row_max); % optimal parameter that yield the best reconstruction result
var2_optimal = var2(col_max);
end

function plotLocal(xy_2compare,grid_GT,metric_GT,sigma,img_h,title_sub)
for j = 1:size(xy_2compare,1)
    xy_fiber = cell2mat(xy_2compare(j));
    plot(xy_fiber(:,2),img_h - xy_fiber(:,1),'Color',[0.5,0.5,0.5],'LineWidth',sigma)
    hold on
end
scatter(grid_GT(:,2),img_h - grid_GT(:,1),sigma^2, metric_GT,'filled')

set(0,'defaultfigurecolor','w')
colormap(flipud(cool))
colorbar
axis square
axis off
title(title_sub)
end

function evaluatePlot(metric_array,x_value,y_value,title_value,dim_img,min_metric,max_metric)
% plot the evaluation results of the two CLAHE parameters under specified kmeans parameter
range_x = [min(x_value) max(x_value)];
range_y = [min(y_value) max(y_value)];

clims = [min_metric max_metric]; % set colormap range
imagesc(range_x,range_y,metric_array,clims)
axis xy % reverse direction of y axis
axis square
% --- set x and y label and ticks -----------------------------------
step_x = (max(x_value) - min(x_value)) / (length(x_value) - 1);
xticks(min(x_value):step_x:max(x_value))
xticklabels(string(x_value))
xlabel('num_c')
ylabel('limit_c')
% --- set title -----------------------------------------------
title(['num_k = ',num2str(title_value,'%.f')])
% ----- set colormap ------------------------------------------
colormap(cool) % set colormap type
colorbar
end
