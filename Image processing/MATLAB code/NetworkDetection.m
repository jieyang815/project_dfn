function [xy_fiber,len_fiber,angle_fiber] = NetworkDetection(img_skel,min_quality)
% --------- Use Depth First Search to trace along fibers ----------------------------------------
% note: each pair of noded is considered to be connected by a fiber
% Input:
% img_skel: skeletonized image
% min_quality: minimum accepted metric value to remove errorneous corners

% Output:
% xy_fiber: location of detected fibers, represents by the coordinate of the pair of endpoints it connects
% len_fiber: array of fiber length
% angle_fiber: array of fiber angle, used to connect broken fibers

xy_fiber = {};
len_fiber = [];
angle_fiber = [];
used_px = false(size(img_skel)); % old pixels to avoid repeating
node_xy = []; % record node coordinates

while sum(sum(used_px)) < sum(sum(img_skel))
    img_unsearched = img_skel - used_px; % pixels not yet searched
    % ---- detect nodes: endpoints and branch points --------------------
    loc_nodes = NodesDetection(img_unsearched,min_quality);
    
    %     % ---- Plot to check ------------------------------------------------------
    %     figure()
    %     imshow(img_unsearched)
    %     hold on
    %     plot(loc_nodes(:,2),loc_nodes(:,1),'gx')
    %     % -------------------------------------------------------------------------
    
    if ~isempty(node_xy)
        node_unused = setdiff(loc_nodes,node_xy,'rows'); % find unused nodes
    else % first round of fiber tracing
        node_unused = loc_nodes;
    end
    node_start = node_unused(1,:); % search starts from the first unused node
    index_node = find((loc_nodes(:,1) == node_start(1))&(loc_nodes(:,2) == node_start(2)));
    queue_px = [node_start node_start index_node]; % queue of pixels to be searched
    
    % ---- Trace fibers along centerlines using DFS ------------------------------
    [xy_fiber,len_fiber,angle_fiber,used_px] = ...
        FILO(queue_px,used_px,img_skel,loc_nodes,len_fiber,angle_fiber,xy_fiber);
    
    node_xy = unique(cell2mat(xy_fiber),'rows'); % node list based on fiber location list
end

%     % ---- plot to check ------------------------------------------------
%     f1 = figure();
%     figure(f1)
%     imshow(img_skel)
%     hold on
%     width_line = 2;
%     for i = 1:length(xy_fiber)
%         temp_fiber = cell2mat(xy_fiber(i));
%         plotFiber(temp_fiber(1,:),temp_fiber(end,:),f1,'y-','bo',width_line)
%     end
%     % ---------------------------------------------------------------------

end

%% Function

function [xy_fiber,len_fiber,angle_fiber,used_px] = ...
    FILO(queue_px,used_px,img_skel,loc_nodes,len_fiber,angle_fiber,xy_fiber)
% Use depth first search to trace fibers, follows the rule: First In Last Out (FILO)

[img_w,img_h] = size(img_skel); % image size & searching boundary
connect_list = []; % connectivity list recording node number on each fiber
while ~isempty(queue_px) % when there are still pixels in queue
    px_current = queue_px(end,1:2); % take out the last pixel in queue
    node_last = queue_px(end,3:4); % last node corresponding to current pixel
    index_last = queue_px(end,5);
    used_px(px_current(1),px_current(2)) = 1; % record used pixels
    queue_px(end,:) = []; % take out the last pixel from queue
    
    px_temp = []; % solid pixels around current pixel
    for i = px_current(1) - 1:px_current(1) + 1 % search 3*3 region around current pixel
        for j = px_current(2) - 1 : px_current(2) + 1
            
            % ----- preliminary check ------------------------------------------------
            if i < 1 || i > img_w || j < 1 || j > img_h % pixel cannot exceed boundary
                continue
            end
            
            % ---- check if current pixel is one of the nodes ---------------------------
            index_current = find((loc_nodes(:,1) == i)&(loc_nodes(:,2) == j),1);
            if ~isempty(index_current)
                if index_current ~= index_last % pixel is one of the nodes, and different from last node
                    px_temp = [px_temp;i j]; % record solid pixels around current search point
                    [node_last,len_fiber,angle_fiber,xy_fiber] = ...
                        createFiber(node_last,[i,j],len_fiber,angle_fiber,xy_fiber);
                    connect_list = [connect_list;index_last index_current]; % add node number to connectivity list
                    index_last = index_current;
                end
                continue
            end
            
            if any((used_px(i,j) == 1))% pixel has been visited
                continue
            end
            
            % ---- Normal solid pixel (on the centerlines but not node) --------
            if img_skel(i,j) == 1 % solid pixel
                px_temp = [px_temp;i j]; % record solid pixel around current pixel
            end
            
        end
    end
    
    [num_row,~] = size(px_temp);
    if num_row ~= 0
        queue_px = [queue_px;px_temp [node_last index_last] .* ones(num_row,3)]; % add solid pixels to searching queue
        
    elseif num_row == 0 && pdist2(px_current, node_last) > 5 % no solid pixel around current pixel & not too close to last node
        dist_temp =  pdist2(px_current, loc_nodes); % calculate its distance to all nodes
        index_list = find(dist_temp < 5);
        
        if  ~isempty(index_list) % close to at least one of the nodes
            % ---- check if connection has been built --------------------
            flag_connect = 0; % initialize flag to mark connection
            if ~isempty(connect_list) % connection has been made before
                for i = 1:length(index_list)
                    index_temp = index_list(i);
                    connect_first = [index_last index_temp];
                    connect_second = [index_temp index_last];
                    if sum(all(connect_list == connect_first,2)) || sum(all(connect_list == connect_second,2)) % conenction has been made
                        flag_connect = 1;
                        break
                    end
                end
            else % no previous connection has been made
            end
            
            if ~flag_connect % no previous or same connection
                [~,index_min] = min(dist_temp); % connect with the closest node
                connect_list = [connect_list;index_last index_min];
                [~,len_fiber,angle_fiber,xy_fiber] = ...
                    createFiber(node_last,loc_nodes(index_min,:),len_fiber,angle_fiber,xy_fiber);
                
            end
            
        else % no nearby detected nodes, then add the current pixel as new node
            [~,len_fiber,angle_fiber,xy_fiber] = ...
                createFiber(node_last,px_current,len_fiber,angle_fiber,xy_fiber);
            loc_nodes = [loc_nodes;px_current];
            connect_list = [connect_list;index_last size(loc_nodes,1)];
            
        end
        %         % ------ Plot current fiber --------------------------------
        %         f0 = figure();
        %         figure(f0)
        %         imshow(img_skel)
        %         hold on
        %         width_line = 2;
        %         plotFiber(loc_nodes(connect_list(end,1),:),loc_nodes(connect_list(end,2),:),f0,'y-','yo',width_line)
        %         % ------------------------------------------------------------
        
    end
end

end

function [node_last,len_list,angle_list,xy_list] = createFiber(node_last,node_current,len_list,angle_list,xy_list)
len_temp = pdist2(node_last,node_current); % calculate fiber length
len_list = [len_list; len_temp]; % add fiber length to global list

slope_temp = (node_last(1) - node_current(1))/(node_last(2) - node_current(2));  % calculate fiber slope
angle_temp = atan(slope_temp) / pi() * 180;
if abs(angle_temp) == 90
    angle_temp = 90;
end
angle_list = [angle_list;angle_temp]; % add fiber slope to global list
xy_list = [xy_list;{[node_last;node_current]}]; % record line location: endpoint coordinates
node_last = node_current; % update last node location
end

%     function plotFiber(node_last,node_current,num_figure,type_line,type_marker,width_line)
%     x_fiber = [node_last(2) node_current(2)];
%     y_fiber = [node_last(1) node_current(1)];
%
%     figure(num_figure)
% %     plot(x_fiber,y_fiber,type_marker) % plot endpoints
% %     hold on
%     plot(x_fiber,y_fiber,type_line,'LineWidth',width_line) % plot fiber
%     hold on
%     end

function loc_nodes = NodesDetection(img_skel,min_quality)
% --- check image type --------------------------------------------------
if ~isa(img_skel,'logical')
    img_skel = logical(img_skel);
end

% ---- Find branch points and endpoints ---------------------------------------
img_branchP = bwmorph(img_skel, 'branchpoints'); % find branch point using morphorlogial operation
img_endP = bwmorph(img_skel, 'endpoints'); % find end point
[row_branchP,col_branchP] = find(img_branchP == 1);
[row_endP,col_endP] = find(img_endP == 1);
xy_eb = [row_endP,col_endP;[row_branchP,col_branchP]]; % xy coordinate of endpoints and branch points

%         % --- Plot nodes on network skeleton to check results ---------------------
%     figure()
%     imshow(img_skel)
%     hold on
%     plot(xy_eb(:,2), xy_eb(:,1),'go','DisplayName','Nodes (morphological)')
%     legend
%     title('Nodes detection')
%     % -------------------------------------------------------------------------

%------- Find corner points (branch points) ------------------------------------------------
corner = myDetectHarrisFeatures(img_skel,'MinQuality',min_quality); % detect crosslinks using Harris corner operator
xy_corner = corner.Location;
xy_corner = uint16(xy_corner);
xy_corner = [xy_corner(:,2) xy_corner(:,1)]; % coordinates of corner points

% ---- merge nodes that are too close -----------------------------------
space = 2;
eb_merged = nodeMerge(xy_eb,xy_eb,space); % merge branchpoints that are close to endpoints
if ~isempty(xy_corner)
    xy_corner = double(xy_corner);
    corner_merged = nodeMerge(xy_corner,xy_corner,space);
    [node_merged,corner_merged2] = nodeMerge(eb_merged,corner_merged,space);
    
    %         % --- Plot nodes on network skeleton to check results ---------------------
    %         figure()
    %         imshow(img_skel)
    %         hold on
    %         %     plot(loc_nodes(:,2), loc_nodes(:,1),'go','DisplayName','Nodes')
    %         eb_merged2 = setdiff(node_merged,corner_merged2,'rows');
    %         plot(eb_merged2(:,2),eb_merged2(:,1),'bx','DisplayName','Nodes (Morphology)')
    %         plot(corner_merged2(:,2), corner_merged2(:,1),'go','DisplayName','Nodes (Harris)','MarkerFaceColor','r')
    %         legend
    %         title('Nodes detection')
    %         % -------------------------------------------------------------------------
    
else
    node_merged = eb_merged;
end

loc_nodes = double(unique(node_merged,'rows')); % convert to double, and remove duplicate nodes

end

function [xy_merged,list_2merge] = nodeMerge(list_2compare,list_2merge,space)
% ---- Merge nodes that are too close -------------------------------------
% Input:
% space: minimum distance between different nodes

% Output
% xy_merged: coordinate of merge nodes
xy_merged = [];

% ---- check if the two node lists are the same -----------------------
if isequal(list_2compare,list_2merge)
    flag_equal = 1;
else
    flag_equal = 0;
end

while ~isempty(list_2compare)
    x_temp =list_2compare(end,2);
    y_temp = list_2compare(end,1);
    xy_merged = [xy_merged;list_2compare(end,:)];
    list_2compare(end,:) = [];
    
    x_polygon = [(x_temp - space) (x_temp + space) (x_temp + space) (x_temp - space)  (x_temp - space)];
    y_polygon = [(y_temp - space) (y_temp - space) (y_temp + space) (y_temp + space)  (y_temp - space)];
    if flag_equal
        index_inside = inpolygon(list_2compare(:,2),list_2compare(:,1),x_polygon,y_polygon); % check if there are close corner points
        list_2compare(index_inside,:) = [];
    else
        index_inside = inpolygon(list_2merge(:,2),list_2merge(:,1),x_polygon,y_polygon); % check if there are close corner points
        list_2merge(index_inside,:) = [];
    end
end

if ~flag_equal
    xy_merged = [xy_merged;list_2merge]; % add the unmerged node into node coordinate list
end

end