function xy_connected = ...
    FiberSeg(img_original,parameter_pre,num_cluster,...
    size_close,len_short,min_quality,error_acc,error_ll)

img_preprocess = imgPreprocess(img_original,parameter_pre);

% %% ---- Binarization-----------------------------------------------
img_binarized = imgBinarize(img_original,img_preprocess,num_cluster,size_close);

% %% ----- Extract skeleton --------------------------------------------
img_skel = imSkel(img_binarized,len_short);

% %% ---- Network detection --------------------------------------------
% ---- use DFS to find fiber ----------------------------------------------
[xy_fiber,len_fiber,angle_fiber] = NetworkDetection(img_skel,min_quality);

% %% Postprocessing: connect broken fibers
%%---- iterate to connect fibers until no new connection is created in one loop %%
max_gap = mean(len_fiber); % tolerance of distance between nearby nodes: average of fiber length
[xy_connected,endP_connected,angle_connected] = ...
    FiberConnect(xy_fiber,angle_fiber,max_gap,img_original,error_acc,error_ll);
end

function img_preprocess = imgPreprocess(img_original,parameter_temp)
% Preprocessing of image
% Input:
% img_original: original image for preprocessing
% parameter_temp: array of preprocessing parameters for CLAHE

% Parameter used in CLAHE:
% limit_clip; contrast limit ([0,1],contrast increase as limit increase ) for contrast-limited adaptive histogram equalization (CLACHE)
% numtiles: image is divided into numtiles*numtiles tiles, numtiles > 2

% Output:
% img_preprocess: image after preprocessing
limit_clip = parameter_temp(1);
numtiles = parameter_temp(2);
img_preprocess = adapthisteq(img_original,'clipLimit',limit_clip,...
    "NumTiles",[numtiles numtiles] );
end

function img_binarized = imgBinarize(img_original,img_preprocess,num_cluster,size_close)
img_binary = segKmeans(img_preprocess,num_cluster);
img_name_bi = 'Binarized (Kmeans)';
% ---- close small holes -----------------
img_binarized = imclose(img_binary,ones(size_close)); % close small holes

% % % ---- Plot binarized image ---------------------------------------------
% % imgCompare(img_original,img_binarized,'Original',img_name_bi)
% img_otsu = imbinarize(img_preprocess);
% imgCompare(img_otsu,img_binarized,'Otsu','Kmeans')
% % % -------------------------------------------------------------------------
end

function img_seg = segKmeans(img_toseg,num_cluster)
%%Use K-means to divide image into k clusters,
%%and extract object from background
[img_cluster,xy_center] = imsegkmeans(img_toseg,num_cluster);

% % ---- Plot to check ---------------------------------------
% figure()
% imshow(labeloverlay(img_toseg,img_cluster))
% title('Kmeans-clustering')
% % ----------------------------------------------------------

[~,cluster_min] = min(xy_center); % background cluster
img_seg = ~(img_cluster == cluster_min);
end

function img_skel = imSkel(img_toskel,len_short)
% --- check image type --------------------------------------------------
if ~isa(img_toskel,'logical')
    img_toskel = logical(img_toskel);
end
% --- Skeletonization --------------------------------------
img_thin = bwmorph(img_toskel,'thin',inf); % Use thinning to skeletonize
img_skel = bwskel(img_thin,'MinBranchLength',len_short); % remove short branches
% imgCompare(img_thin ,img_skel,'Thinning','Small branches removed')

% % ---- plot skeletonized image -------------------------------------------
% imgCompare(img_toskel ,img_skel,'Binarized','Skeletonized')
skel_parallel = bwskel(img_toskel,'MinBranchLength',len_short);
% imgCompare(skel_parallel,img_skel,'Parallel','Sequantial')
% % -------------------------------------------------------------------------

% f8 = figure();
% figure(f8)
% % imshow(imresize(img_skel,10))
% imshow(img_skel)
end

function imgCompare(img1,img2,title_img1,title_img2)
figure()
subplot(1,2,1)
imshow(img1,[])
title(title_img1)
subplot(1,2,2)
imshow(img2)
title(title_img2)
end

