%--------------------------------------------------------------------------
% Metric to evaluate the matching rate between ground truth and segmentaion result
% -------------------------------------------------------------------------
function [Pmean,Nmean,Pmetric_GT,Nmetric_SEG,grid_GT,grid_seg] = ...
    getMetric(xy_GT,xy_seg,sigma,spacing)
% Input:
% xy_GT: coordinate of fibers as ground truth
% xy_seg: coordinate of fibers to evaluate
% sigma: sensitive parameter
% spacing: sampling internal

% Output:
% Pmean: average of the positive metric
% Nmean: average of the negative metric
% Pmetric_GT: positive gaussian metric of ground truth, defined as the rate between
%                length of correctly detected fibers and length of 
%                all fibers in ground truth
% Nmetric_SEG: negative gaussian metric of segmentation result, defined as the rate between
%                length of falsely detected fibers and length of 
%                all fibers in segmentation results

% ---- Sample along the centerline of the fibers on a uniform grid --------------
grid_GT = uniformSample(xy_GT,spacing);
grid_seg = uniformSample(xy_seg,spacing);

% ----- Calculate gaussian metrics between two networks -------------------
type_metric = 1; % 1: positive mean metric; 2: negative mean metric
[Pmetric_GT,Pmean] = gaussianMetric(grid_seg,grid_GT,sigma,type_metric); % Fibers exist in ground truth, but missed in segmentation
type_metric = 2; % 1: positive mean metric; 2: negative mean metric
[Nmetric_SEG,Nmean] = gaussianMetric(grid_GT,grid_seg,sigma,type_metric); % Fibers detected in segmentation, but did not exist in ground truth
end

%% Functions
function xy_grid = uniformSample(xy_data,spacing)
% ---- Sample the fiber with fixed interval ------------------------------
% Input:
% xy_data [cell array]: fiber coordinates, consists of 2 or more nodes
% spacing: sampling interval, which is the grid size

% Output:
% xy_grid: coordinates of points on the grid
xy_grid = [];
for i = 1:size(xy_data,1)
    xy_temp = cell2mat(xy_data(i));
    for j = 1:size(xy_temp)-1
        xy_start = xy_temp(j,:);
        xy_end = xy_temp(j + 1,:);
        len_temp = pdist([xy_end(2) xy_end(1);
                          xy_start(2) xy_start(1)]);
        num_spacing = ceil(len_temp / spacing); % calculate number of intervals
        sampling_vector = (0:1/num_spacing:1)';
        x_sampling = xy_start(2) + (xy_end(2) - xy_start(2)) .* sampling_vector;
        y_sampling = xy_start(1) + (xy_end(1) - xy_start(1)) .* sampling_vector;
        xy_grid = [xy_grid;y_sampling x_sampling];
    end
end
end

function [gaussian_X2Y,mean_metric] = gaussianMetric(X,Y,sigma,type_metric)
% -------------------------------------------------------------------------
% Finds the nearest neighbor in X for each query point in Y  
% and calculated distances between query point and its neighbor
% finally the distances are converted to gaussian metric
% ------------------------------------------------------------------------
% Input:
% sigma: sensitivity parameter for gaussian
% X: data point to evaluate
% Y: ground truth

% Output:
% gaussian_X2Y: gaussian metric array
% mean_metric: average of gaussian metric

index_neighbor = knnsearch(X,Y); % find the index of nearest neighbor in X for query point in Y
xy_neighbor = X(index_neighbor,:); % corresponding neighbor coordinates

gaussian_X2Y = zeros(size(Y,1),1);
for i = 1:size(Y,1)
    xy_2query = Y(i,:);
    xy_ref = xy_neighbor(i,:);
    dist_temp = pdist([xy_2query; xy_ref]);
    gaussian_X2Y(i) = exp(-0.5 * ((dist_temp^2) / (sigma^2)) ); % convert distance into gaussian metric
end

switch type_metric
    case 1 % positve metric
        mean_metric = mean(gaussian_X2Y);
    case 2 % negative metric
        mean_metric = 1 - mean(gaussian_X2Y);
end

end