function [node_matrix,connect_elem,n2_elem,node_BC,disp_BC,coord_BC] = createBeamMesh(geo_fiber,len_elem,tol_distBC,strain)
% Create mesh for the beam element
% ---------------------------------------------------------------
% Input:
% ax_beam: beam axis, represents by two endpoints
% order_beam: type of beam

% Output:
% node_matrix: index and coordinates of the element node [index1 x1 y1 z1;index2 x2 y2 z2]
% connect_element: connectivity information of element [index_elem1,indexEP1_eleml indexEP2_elem1;
%                                                       index_elem2,indexEP21_elem2 indexEP2_elem2]
% n2_elem: normal vector of element [x_v1 y_v1 0;x_v2 y_v2 0]
% node_BC: boundary node set, record the index of boundary (struct - Xmin, Xmax,Ymin,Ymax)
% coord_BC: record the coordinate of nodes at the boundary
% ---------------------------------------------------------------
ax_beam = geo_fiber.ax_beam;
len_beam = geo_fiber.len_beam;
n2_beam = geo_fiber.n2_beam;
node_xyz = [];
elem_xyz = {};
n2_elem = [];
for i = 1:length(ax_beam)
    % --- import geometry information of current beam element -------------
    temp_beam = ax_beam(i,:); 
    temp_n2 = n2_beam(i,:);
    temp_len = len_beam(i); 
    
    xyz_ep1 = cell2mat(temp_beam(1)); % coordinate of first endpoint
    xyz_ep2 = cell2mat(temp_beam(2));

    temp_num = ceil(temp_len / len_elem); % calculate number of intervals
    vector_mesh = (0:1/temp_num:1)'; % mesh vector/sampling vector
    x_elem = xyz_ep1(1) + (xyz_ep2(1) - xyz_ep1(1)) .* vector_mesh;
    y_elem = xyz_ep1(2) + (xyz_ep2(2) - xyz_ep1(2)) .* vector_mesh;
    
%     % ----------- plot to check discretization --------------------------------------
%     plot([xyz_ep1(1) xyz_ep2(1)],[xyz_ep1(2) xyz_ep2(2)],'Color',[0.5,0.5,0.5],'LineWidth',5)
%     hold on
%     plot(x_elem,y_elem,'ro','MarkerFaceColor','r')
%     % ----------------------------------------------------------------
    
    % --- record element coordinate for constructing connectivity matrix ---
    count_elem = 0; % count number of element in current beam
    for j = 1:(length(x_elem) - 1)
        elem_ep1 =  [x_elem(j) y_elem(j) 0];
        elem_ep2 =  [x_elem(j + 1) y_elem(j + 1) 0];
        elem_xyz = [elem_xyz;{elem_ep1 elem_ep2}];
        count_elem = count_elem + 1;
    end
    % ----- record element normal vector ------------------------
    n2_elem = [n2_elem;repmat(temp_n2,count_elem,1)]; 
    % --- update node coordinate --------------------------------
    node_xyz = [node_xyz;x_elem y_elem zeros(size(x_elem))]; % record node coordinate
end
% --- arrange node matrix and remove duplicate node -----------------------
node_xyz = unique(node_xyz,'rows');
node_index = (1:size(node_xyz,1))';
node_matrix = [node_index node_xyz];

% ----- set up connectivity matrix -------------------------------------
connect_elem = zeros(size(elem_xyz,1),3);
for i = 1:size(elem_xyz,1)
   elem_ep1 = cell2mat(elem_xyz(i,1));
   elem_ep2 = cell2mat(elem_xyz(i,2));
   [~,idx_ep1] = ismember(node_xyz,elem_ep1,'rows'); % find the index of node
   [~,idx_ep2] = ismember(node_xyz,elem_ep2,'rows');
    connect_elem(i,:) = [i find(idx_ep1 == 1) find(idx_ep2 == 1) ]; % add index of the two endpoints to connectivity matrix          
end

% -------- set up node set at the boundary -------------------------------
lim_x = max(node_xyz(:,1)); % x and y limit for plotting
lim_y = max(node_xyz(:,2));
idx0_Xmin = find(node_xyz(:,1) <= min(node_xyz(:,1)) + tol_distBC); % left boundary
idx0_Xmax = find(node_xyz(:,1) >= max(node_xyz(:,1)) - tol_distBC);
idx_Ymin = find(node_xyz(:,2) <= min(node_xyz(:,2)) + tol_distBC);
idx_Ymax = find(node_xyz(:,2) >= max(node_xyz(:,2)) - tol_distBC);

% --- remove duplicate node -----------------------------------------------
idx1_Xmin = setdiff(idx0_Xmin,idx_Ymin); 
idx_Xmin = setdiff(idx1_Xmin,idx_Ymax); % assign corner node to top boundary
idx1_Xmax = setdiff(idx0_Xmax,idx_Ymin);
idx_Xmax = setdiff(idx1_Xmax,idx_Ymax); % assign corner node to bottom boundary

% ---- save index of boundary node set ------------------------------------
node_BC.Xmin = idx_Xmin;
node_BC.Xmax = idx_Xmax;
node_BC.Ymin = idx_Ymin;
node_BC.Ymax = idx_Ymax;

% ---- save coordinate of boundary node set -------------------------------
xy_left = node_xyz(idx_Xmin,1:2);
xy_right = node_xyz(idx_Xmax,1:2);
xy_lower = node_xyz(idx_Ymin,1:2);
xy_upper = node_xyz(idx_Ymax,1:2);
coord_BC.left = [idx_Xmin xy_left];
coord_BC.right = [idx_Xmax xy_right];
coord_BC.lower = [idx_Ymin xy_lower];
coord_BC.upper = [idx_Ymax xy_upper];

% % ---- Plot to check -------------------------------------------------------
% plot(xy_left(:,1),xy_left(:,2),'go','MarkerFaceColor','g','MarkerSize',5)
% plot(xy_right(:,1),xy_right(:,2),'go','MarkerFaceColor','g','MarkerSize',5)
% plot(xy_lower(:,1),xy_lower(:,2),'bo','MarkerFaceColor','b','MarkerSize',5)
% plot(xy_upper(:,1),xy_upper(:,2),'bo','MarkerFaceColor','b','MarkerSize',5)
% xlim([0 lim_x + 1])
% ylim([0 lim_y + 1])
% axis off
% % --------------------------------------------------------------------------

% ------ defined boundary condition ---------------------------------------
width = max(node_xyz(:,1)) - min(node_xyz(:,1)); % size of the fiber network
disp_BC.Xleft = - (strain * width) / 2; %Displacement at left boundary
disp_BC.Xright = (strain * width) / 2; %Displacement at right boundary
disp_BC.Yupper = (strain * width) / 2; %isplacement at upper boundary
disp_BC.Ylower = - (strain * width) / 2; %Displacement in lower boundary

end