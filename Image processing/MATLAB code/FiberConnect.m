function [xy_connected,endP_connected,angle_connected] = ...
    FiberConnect(xy_fiber,angle_fiber,max_gap,~,error_acc,error_ll)
% ------- Find connection between broken fibers -----------------------------------------
% Input:
% xy_list: fiber location
% anglie_list: fiber angle
% max_gap: distance criteria to check if a fiber is nearby
% error_acc: angle criteria to check if two fibers are of similar angle
% error_ll: angle criteria to check if two fibers is parallel

% Output:
% xy_connected: fiber location, including all nodes on a fiber
% endP_connected: endpoints of connected fibers
% angle_connected: angle of connected fibers

% ------ preliminary setup ----------------------------------------------------------------------------
xy_unconnect = xy_fiber;
endP_unconnect = xy_fiber; % cell array: endpoints of fibers
angle_unconnect = angle_fiber;
flag_isnew = 1; % flag to mark if new connection is made in major loop
flag_continue = 0; % flag to mark if new connection has been made in minor loop

while flag_isnew % new connection is made in last major loop
    % ---- Preliminary setup ------------------------------------------------------
    xy_2connect = xy_unconnect; % coordinate list for fibers to be connected in current major loop
    angle_2connect = angle_unconnect;
    endP_2connect = endP_unconnect;
    xy_unconnect = []; % new cell arrary for fibers to be connected in the next major loop
    angle_unconnect = [];
    endP_unconnect = [];
    flag_isnew = 0;
    
    % ----- minor loop to connect fibers -----------------------------
    while ~isempty(endP_2connect)
        if ~flag_continue % no new connection from last minor loop
            xy_temp = cell2mat(xy_2connect(end)); % take one fiber from unconnect list
            angle_temp = angle_2connect(end); % angle of current fiber
            endP_temp = cell2mat(endP_2connect(end));
            xy_2connect(end) = []; % remove current fiber from unconnect list
            endP_2connect(end) = [];
            angle_2connect(end) = [];
        end
        
        flag_continue = 0;
        
        % ---- 1. Select fibers with similar slope ------------------------------- %
        min_angle = angle_temp - error_acc; % angle range
        max_angle = angle_temp + error_acc;
        logical_angle = (angle_2connect > min_angle) & (angle_2connect < max_angle);
        if ~any(logical_angle) % no fibers with similar slope exist ---> no new connection to be made
            
            xy_unconnect = [xy_unconnect;{xy_temp}]; % add current fiber to unconnected list
            angle_unconnect = [angle_unconnect;angle_temp];
            endP_unconnect = [endP_unconnect;{endP_temp}];
            
            continue
        else % fibers with similar orientation exist
            xy_angle = xy_2connect(logical_angle); % location of similar fibers
            endP_angle = endP_2connect(logical_angle);
            index_angle = find(logical_angle == 1); % fiber index in list_unconnect, for later removal
        end
        
        % ----- 2. Select fibers that are close -----------------------------------------------------------
        dist_node = (pdist2(endP_temp,cell2mat(endP_angle)))'; % calculate distance between current fiber and other fibers in list
        logical_close = dist_node < max_gap; % check if any fiber with similar slope is also close
        if sum(sum(logical_close)) == 0 % no close fibers
            
            xy_unconnect = [xy_unconnect;{xy_temp}];
            angle_unconnect = [angle_unconnect;angle_temp];
            endP_unconnect = [endP_unconnect;{endP_temp}];

            continue
        end
        
        % ----- set up node list of nearby fibers --------------------------
        [endP_close,endP_far,index_2keep,flag_col,xy_close,dist_sorted] = ...
            FindCloseNode(logical_close,xy_angle,endP_angle,index_angle,dist_node);
        
        % ---- 3. Check if close fibers are parallel to current fibers ------------------
        flag_2remove = zeros(size(index_2keep,1),1); % record logical index of connected fiber
        for i = 1:size(endP_close,1)
            
            endP_parallel = [endP_close(i,:);endP_far(i,:)];
            xy_parallel = xy_close(i);
            
            % set up angle tolerance for parallel condition
            if dist_sorted(i) < max_gap / 2
                error_check = error_acc;
            elseif dist_sorted(i) >= max_gap / 2 % fibers that are not so close should have lower angle tolerance
                error_check = error_ll;
            end
            min_angle = angle_temp - error_check;
            max_angle = angle_temp + error_check;
            % ------- check if parallel condition is fulfilled ------------------------
            [flag_parallel,~] = isParallel(min_angle,max_angle,endP_temp,endP_parallel);            
            if flag_parallel % parallel
                continue
            end
            
            % --- if fiber is not parallel, then connect it with current fiber ----------------------
            flag_2remove(i) = 1; % record location of connected fibers in index_2check
            flag_continue = 1; % flag to mark if new fiber is create in minor loop
            flag_isnew = 1; % new connection is made in major loop
            xy_temp = [xy_temp;cell2mat(xy_parallel)]; % add node coordinate of the nearby fibers to create new fiber
            xy_temp = unique(xy_temp,'rows'); % avoid duplicate nodes
            
            % ---- arrange order of nodes on the fiber -----------------------------
            if flag_col(i) == 1 % the first endpoint on current fiber is connected
                % the further endpoint on another fiber becomes the new first endpoint
                dist_2first = (pdist2(endP_parallel(2,:),xy_temp));
                [~,index_sorted] = sort(dist_2first,'ascend'); % sort according to node distance to the first endpoint
            elseif flag_col(i) == 2 % the second endpoint on current fiber is connected
                dist_2last = (pdist2(endP_parallel(2,:),xy_temp));
                [~,index_sorted] = sort(dist_2last,'descend'); % sort according to node distance to the second endpoint
            end
            xy_temp = xy_temp(index_sorted,:); % sort coordinate of nodes according to distance to the connected endpoint
            endP_temp = [xy_temp(1,:);xy_temp(end,:)];
            slope_temp = (endP_temp(1,1) - endP_temp(2,1)) /  (endP_temp(1,2) - endP_temp(2,2));
            angle_temp = atan(slope_temp) / pi() * 180; % update fiber angle
            
        end
        
        % --- Postprocessing after filtering and connecting fibers -----------------------------
        
        % ----- 1. remove connected fiber from unconnect list -----------------------
        if flag_continue
            index_2remove = index_2keep(logical(flag_2remove));
            endP_2connect(index_2remove,:) = [];
            xy_2connect(index_2remove,:) = []; % coordinate list for fibers to be connected in current loop
            angle_2connect(index_2remove,:) = [];
        end
        
        % ---- 2. check if no new connection is created or last fiber has been taken out from list
        if ~flag_continue || isempty(endP_2connect)
            xy_unconnect = [xy_unconnect;{xy_temp}];
            angle_unconnect = [angle_unconnect;angle_temp];
            endP_unconnect = [endP_unconnect;{endP_temp}];           
            continue           
        end
        
    end
    
end

% ----- Set up connected fiber list ---------------------------------------
xy_connected = xy_unconnect;
endP_connected = endP_unconnect;
angle_connected = angle_unconnect;

end


%% Function
function [flag_parallel,angle_2check] = isParallel(min_angle,max_angle,xy_current,xy_2connect)
% --- function to check if fiber is parallel ------------------------------
% Select one node from two fibers separately, and calculate to see if the
% line between the pair of the node is of the same slope as current fiber
% if the slope of every line is the within the angle tolerance, then it is not parallel to current fiber
count_inrange = 0;
count = 0;
for i = 1:2
    for j = 1:2
        x1 = xy_current(i,2);
        y1 = xy_current(i,1);
        x2 = xy_2connect(j,2);
        y2 = xy_2connect(j,1);
        if x1 == x2 && y1 == y2 % share the same endpoint
            continue
        end
        slope_2check = (y2 - y1) / (x2 - x1);
        angle_2check = atan(slope_2check) / pi() * 180;
        count = count + 1;
        if angle_2check > min_angle && angle_2check < max_angle % within tolerance
            count_inrange = count_inrange + 1;
        end
    end
end
if count == count_inrange % compare number of lines that fulfill the angle tolerance to total number
    flag_parallel = 0;
else
    flag_parallel = 1;
end
end

function [endP_close,endP_far,index_2keep,flag_col,xy_close,dist_sorted] = ...
          FindCloseNode(logical_close,xy_angle,endP_angle,index_angle,dist_node)

% --- Check which pair of nodes are close, and set up a node list containing nodes that are closer
% and another list with nodes that are further -------------

% Input:
% logical_close: logical index of fibers in endP_angle that are nearby
% xy_angle: coordinate of fibers that are of similar angle as current fiber
% endP_angle: endpoints of fibers that are of similar angle as current fiber
% index_angle: index of close fibers in unconnect fier list
% dist_node: distance between nodes

% Output:
% endP_close: close nodes on the nearby fibers
% endP_far: far nodes on the nearby fibers
% index_2keep: index of fibers that are kept for later fiber filtering
% flag_col: index of the node on current fiber that is close to nodes on other fibers
% xy_close [cell array]: coordinate of close fibers
% dist_sorted: sorted distance between close nodes

index_close = find((logical_close(:,1) + logical_close(:,2)) > 0); % index of close node
endP_close = []; % endpoints that are closer to current fiber
endP_far = []; % endpoints that are further to current fiber
index_temp = 0;
index_xyangle = []; % index of nearby fibers (/fibers to keep) in index_toremove
dist_2sort = []; % distance between a pair of close endpoints, the fiber connection will start from the closest
flag_col = []; % mark index of endpoint to be connected at current fiber
for i = 1:size(index_close,1)
    % --- check if a pair of nodes is used --------------------------
    if i ~= 1 && ((index_temp + 1) == index_close(i))
        continue
    end
    
    index_temp = index_close(i); % fiber index in xy_angle
    
    if rem(index_temp,2) ~= 0 % first node at the fiber
        index_xyangle = [index_xyangle;(index_temp + 1) / 2]; % index of close fibers in xy_angle
        endP_2check = cell2mat(endP_angle(index_xyangle(end)));
        % --- check if both nodes are close ---------------------
        if i == size(index_close,1) || index_close(i + 1) ~= index_temp + 1
            % last element in the array or only one node is close to current fiber
            endP_close = [endP_close;endP_2check(1,:)];
            endP_far = [endP_far;endP_2check(2,:)];
            
            [dist_min,col_min] = min(dist_node(index_temp,:));
            dist_2sort = [dist_2sort;dist_min]; % record distance between fibers for sorting
            flag_col = [flag_col;col_min]; % mark endpoint on current fiber
        elseif  index_close(i + 1) == index_temp + 1 % another node is also close to current fiber
            dist_2compare = dist_node(index_temp:index_temp + 1,:);
            min_2compare = min(min(dist_2compare));
            [row_min,col_min] = find(dist_2compare == min_2compare,1);
            % for parallel condition, there might be two same min distance
            % it will be excluded later
            flag_col = [flag_col;col_min];
            
            if row_min == 1 % first node is the closest
                endP_close = [endP_close;endP_2check(1,:)];
                endP_far = [endP_far;endP_2check(2,:)];
            elseif row_min == 2
                endP_close = [endP_close;endP_2check(2,:)];
                endP_far = [endP_far;endP_2check(1,:)];
            end
            dist_2sort = [dist_2sort;min_2compare];
        end
        
    else % second node in a fiber
        index_xyangle = [index_xyangle;index_temp / 2];
        endP_2check = cell2mat(endP_angle(index_xyangle(end)));
        % second node is closest
        endP_close = [endP_close;endP_2check(2,:)];
        endP_far = [endP_far;endP_2check(1,:)];
        
        % min distance and col
        [dist_min,col_min] = min(dist_node(index_temp,:));
        dist_2sort = [dist_2sort;dist_min];
        flag_col = [flag_col;col_min];
    end
end

index_2keep = index_angle(index_xyangle);
xy_close = xy_angle(index_xyangle);

% sort fiber coordinate in ascending order according to their distance to current fiber
if size(dist_2sort,1) > 1
    [dist_sorted,sorted] = sort(dist_2sort,'ascend');
    endP_close = endP_close(sorted,:);
    endP_far = endP_far(sorted,:);
    index_2keep =  index_2keep(sorted);
    flag_col = flag_col(sorted);
    xy_close = xy_close(sorted);
else
    dist_sorted = dist_2sort;
end
end

