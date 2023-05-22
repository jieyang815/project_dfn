clc
clear
%% Preprocess of geometry
% 1. divide a fiber into multiple fiber segments
% 2. remove isolated fibers
% ---------------------------------------------------------------------------

% --- Import geometry of reconstructed fiber network ----------------------
path_data = uigetdir(' ','Select reconstructed fiber network data folder'); % select the data folder
pattern_data = fullfile(path_data,'*.mat');
path_original = dir(pattern_data); % get list of data path

% --- select folder to save preprocessed geomtry --------------------------
path_tosave = uigetdir(' ','Select the folder to save figures');

% -----------------------------------
tic
% -----------------------------------
for i = 1:length(path_original)
    fprintf('Preprocessing for DFN %i begins\n',i)
    
    name_temp = path_original(i).name;
    path_temp = fullfile(path_data, name_temp); % construct complete file path
    data_temp = load(path_temp); % import fiber network data
    xy_list = data_temp.xy_SEG;
    
%     % ---- plot to check -------------------------
%     figure()
%     % ---------------------------------------------

    % ---- record the nodes on the fiber ----------------------------------
    xy_segment = [];
    list_node = [];
    for j = 1:size(xy_list,1)
        xy_fiber = cell2mat(xy_list(j));
        list_node = [list_node;xy_fiber(:,2) xy_fiber(:,1)];
        
%         % --- Plot to check -----------------------------------
%         plot(xy_fiber(:,2),98 - xy_fiber(:,1),'-','Color',[0.5,0.5,0.5],...
%             'LineWidth',3)
%         hold on
%         plot([xy_fiber(1,2) xy_fiber(end,2)],98 - [xy_fiber(1,1) xy_fiber(end,1)],'go')
%         % ---------------------------------------------------
        
        % --- iterate to divide the fiber to multiple segments
        for k = 1:size(xy_fiber,1) - 1
            node_1 = [xy_fiber(k,2) xy_fiber(k,1)];
            node_2 = [xy_fiber(k + 1,2) xy_fiber(k + 1,1)];
            xy_segment = [xy_segment;node_1 node_2]; % record the two endpoints of the segment
        end
    end
    
    % --- construct fiber segment connectivitiy list to remove isolated fibers ---------------------------
    unique_node = unique(list_node,'rows'); % set up node list
    idx_connect = zeros(size(xy_segment,1),2);
    for k = 1:size(xy_segment,1)
        temp_endP1 = xy_segment(k,1:2);
        temp_endP2 = xy_segment(k,3:4);
        [~,idx_ep1] = ismember(unique_node,temp_endP1,'rows'); % find the index of node
        [~,idx_ep2] = ismember(unique_node,temp_endP2,'rows');
        idx_connect(k,:) = [find(idx_ep1 == 1) find(idx_ep2 == 1) ]; % record index of the two endpoints
    end
    temp_edges = 1:size(unique_node,1); % node number
    temp_count = histc(idx_connect(:), temp_edges);
    index_single = find(temp_count == 1); % endpoints appear only once
    logical_single = zeros(size(idx_connect,1),1); % find the index of isolated fibers
    for k = 1:size(index_single,1)
        temp_single = index_single(k);
        [row_temp,~] = find(idx_connect == temp_single);
        logical_single(row_temp) = logical_single(row_temp) + 1;
    end
    idElem_single = find(logical_single > 1); % remove fiber segment with two single endpoints
    
%     % --- plot single fiber to check --------------------------------
%     xy_single = xy_segment(idElem_single,:);
%     for k = 1:size(xy_single,1)
%         temp_endP1 = xy_single(k,1:2);
%         temp_endP2 = xy_single(k,3:4);
%         plot([temp_endP1(1) temp_endP2(1)],98 - [temp_endP1(2) temp_endP2(2)],'r-','LineWidth',3)
%     end
%     axis off
%     % ----------------------------------------------------------------

    xy_segment(idElem_single,:) = []; % remove isolated fibers
    save_path = fullfile(path_tosave,['Geometry_DFN_img',num2str(i,'%.f'),'.mat']);
    save(save_path,'xy_segment')
    
    fprintf('Geometry for Image %i saved\n',i)
end
% --------------
toc
% --------------