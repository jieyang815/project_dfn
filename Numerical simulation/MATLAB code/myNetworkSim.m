% Main Script for importing reconstructed fiber network to Abaqus
% Input:
%     1) Geometry of reconstructed fiber network
%     2) Discretization Parameters
%     3) Material Properties
%
% Output:

%     4) Abaqus input file

%% INPUT
clc
clear
% --- Predefined ----------------------------------------------------------
% ---- Beam Parameters ----------------------------------------------------
order_beam = 1; % B21:timensinko beam
dim_beam = 2; % dimension of beam
r_beam = 1; % radius of the beam (micrometer)
youngs_beam = 100; % 100 Mpa = 100 uN/um2
nu_beam = 0.45; % Poisson's Ratio for beams
param_beam = [youngs_beam nu_beam];
n1_elem = [0,0,-1]; % normal vector of 2d element
scale_px2um = 2; % scaling factor from pixel to micrometer

% --- Discretization: quasi distant element -------------------------------
len_elem = 1;

% ---- Boundary node set --------------------------------------------------
num_BC = 3; % 
tol_distBC = num_BC * len_elem; % node lie within this distance will be considered as boundary node set

% ----- Dsiplacement Mode -------------------------------------------------
disp_mode = 'BiaxialExt'; % Biaxial extension
strain = 0.125; % strain 

% Abaqus Steps
nSteps = 20;
param_step. InitStep = 0.01;
param_step. minStep = 1e-5;
param_step. maxStep = 1;
param_step. STBLfac = 2e-4;

%% ------- Load Fiber Data (Straight fiber segments) ------------------------
path_data = uigetdir(' ','Select fiber segment folder'); % select the data folder
pattern_data = fullfile(path_data,'*.mat');
path_original = dir(pattern_data); % get list of data path

%% ---- Select directory to save input file ---------------------------
dir_inp = uigetdir(' ','Select folder to save Abaqus input file'); % select the data folder

%% Discretization and implementation of abaqus input file
list_name = [];
list_jobname = [];
for i = 1:length(path_original)
    
    % --- 1. Import fiber segment geometry ---------------------------------
    name_temp = path_original(i).name;
    path_temp = fullfile(path_data, name_temp); % construct complete file path
    name_split = split(char(name_temp),["_","."]); % cut before file number for file name reconstruction
    name_data = string(name_split(2)); % construct file name
    data_Fmodel = load(path_temp); % load fiber geometry data
    segment_xy = double(data_Fmodel.xy_segment)/scale_px2um; % import fiber segment loacation data
    
    % ----- 2. Define beam element -----------------------------------------------
    ax_beam = {}; % beam axis vector
    n2_beam = []; % local 2-axis
    len_beam = []; % length
    ori_beam = []; % orientation
    for j = 1:length(segment_xy)
        if dim_beam == 2 % 2d planar beam
            beam_ep1 = [segment_xy(j,1:2) 0]; % coodinate of node in xy-plane: [x y 0]
            beam_ep2 = [segment_xy(j,3:4) 0];
        elseif dim_beam == 3 % 3d beam
            beam_ep1 = segment_xy(j,1:3);
            beam_ep2 = segment_xy(j,4:6);
        end
        
        % --- geometry of beam -------------
        ax_beam = [ax_beam;{beam_ep1 beam_ep2}]; % beam axis: n*2, 2 endpoints in different column
        len_temp = pdist2(beam_ep1,beam_ep2); % calculate length of fiber segment
        len_beam = [len_beam;len_temp];
        angle_temp = atan((beam_ep2(2) - beam_ep1(2))/...
                    (beam_ep2(1) - beam_ep1(1))) * 180 / pi(); % ---- calculate orientation of the fiber
        if angle_temp < 0 % > 90
            angle_temp = 180 - abs(angle_temp);
        end
        ori_beam = [ori_beam;angle_temp];
        tangent_temp = beam_ep2 - beam_ep1; % tangent vector of the cross section
        n2_temp = cross(tangent_temp,n1_elem); % 2-axis = cross product of tangent vector and n1
        n2_beam = [n2_beam;n2_temp];
    end
    % ----- add to geometry data set --------------------------------------
    geo_fiber.ax_beam = ax_beam;
    geo_fiber.len_beam = len_beam;
    geo_fiber.n2_beam = n2_beam;    
%     save(['Geo_v4_',char(name_data),'.mat'],'geo_fiber')
    
    fprintf('Geometry import completed for image %i. Average segment length is %.2f, mean orientation is %.2f\n',...
        i,mean(len_beam))
    
     
%     % ------- Plot to check the discretization result or boundary node set 
%     f1 = figure();
%     figure(f1)
%     % --------------------------------------------------------------------
        
    % ---- Discretization of geometry ------------------------------------
    fprintf('Generating BEAM mesh...\n')       
    [node_matrix,connect_elem,n2_elem,node_BC,disp_BC,coord_BC] = ...
        createBeamMesh(geo_fiber,len_elem,tol_distBC,strain);

    % --- save the dicretization result and coordinate of boundary node set (for postprocessing) ----------------------------------
%     save([dir_inp,'\Mesh_se012_',char(name_data),'.mat'],'node_matrix','connect_elem','n2_elem','node_BC')
    save([dir_inp,'\Coord_initial_',char(name_data),'.mat'],'coord_BC')

    fprintf('Mesh generation completed for image %i \n',i)
    
    % -------- writing abaqus input file --------------------------------------
    fprintf('Writing Abaqus input file for image %i ...\n',i)
    filename_inp = writeINP(name_data,dir_inp,node_matrix,node_BC,n1_elem,n2_elem,connect_elem,...
                            order_beam,i,r_beam,param_beam,param_step,disp_BC);

    fprintf('\n Abaqus input file for image %i completed ...\n',i)
end
%% Function
function inputCheck(data_test,file_inp,num_item)
% - check if length of the array exceed the limit for each line in inp file --
len_temp = length(data_test); % length of the array
if len_temp < num_item
    str_line = ['%i', repmat(', %i',1,len_temp - 1), ', \n'];
    fprintf(file_inp, str_line, data_test);
elseif len_temp > num_item
    num_row = floor(len_temp/num_item);
    str_line = ['%i', repmat(', %i',1,num_item - 1), ', \n'];
    fprintf(file_inp, str_line, data_test(1:num_row*num_item));
    % --- check if length of the data is the multiple of num_item -------
    rem_num = rem(len_temp,num_item); % calculate the remainder
    if rem_num > 0
        str_line = ['%i', repmat(', %i',1,rem_num - 1), '\n']; % write the remaining lines
        fprintf(file_inp, str_line, data_test(end - rem_num + 1:end));
    end
end
end

function filename_inp = writeINP(name_data,dir_inp,node_matrix,node_BC,n1_elem,n2_elem,connect_elem,...
    order_beam,i,r_beam,param_beam,param_step,disp_BC)
% Create Abaqus input file ------------------------------------------------
% Input:
% beam geometry
% boundary condition
% material parameter

% Output
% Abaqus input file
% ------------------------------------------------------------------------

filename_inp = sprintf('%s%s%s','FiberModel_se012unit_',name_data,'.inp');
fullpath_inp = fullfile(dir_inp, filename_inp); % construct complete file path
file_inp = fopen(fullpath_inp,'w'); % create the input abaqus file

% Write NODES matrix
fprintf(file_inp, '*NODE \n'); % write the node coordinate, [NodeNumber x y]
fprintf(file_inp,'%i, %2.14e, %2.14e, %2.14e \n', node_matrix');
% \t: align the output of the node coordinates
% middle part: indicates that each row of data will contain four values, separated by commas
% %i is for an integer value (the node number),
% while %2.14e is for a floating-point value (the x, y, or z coordinate)
% The node_matrix variable contains the nodal coordinates for each node,
% with each row of the matrix representing one node.

% ---- Write Node set at the boundary --------
num_item = 8; % number of item per line (maximum is 16, but comma also counts)
% ---- right boundary
fprintf(file_inp, '*NSET, NSET = leftNSET \n');
inputCheck(node_BC.Xmin,file_inp,num_item) % check the length of array before writing into the file

% ---- left boundary
fprintf(file_inp, '*NSET, NSET = rightNSET \n');
inputCheck(node_BC.Xmax,file_inp,num_item)

% ----- upper boundary
fprintf(file_inp, '*NSET, NSET = lowerNSET \n');
inputCheck(node_BC.Ymin,file_inp,num_item)

% ------ lower boundary
fprintf(file_inp, '*NSET, NSET = upperNSET \n');
inputCheck(node_BC.Ymax,file_inp,num_item)

% ----- Define elements --------------------------------------------------
fprintf(file_inp, '*ELEMENT, TYPE=B2%d, ELSET=BeamElem%d \n', [order_beam,i]);
% ELSET is the name of the element set, B2%d is replaced by order, and last %d is replaced by i
if order_beam == 1 || order_beam == 3 % 2 node beam element
    fprintf(file_inp, ['%i', repmat(', %i',1,2), ' \n'], connect_elem');
    %  element number and two node numbers
else
    fprintf(file_inp, ['%i', repmat(', %i',1,3), ' \n'],connect_elem');
end

%Define Normals
fprintf(file_inp, '*NORMAL,TYPE=ELEMENT\n');
fprintf(file_inp, '%d, %d, %.6f, %.6f, %.6f \n', [connect_elem(:,1:2),n2_elem]');
% the normals of the beam elements: %d represents the element node numbers, ...
% and %f represents the normal direction values: x,y,z

% Define Cross sections
fprintf(file_inp, '*BEAM SECTION, ELSET=BeamElem%d, MATERIAL=BeamMatl, SECTION=CIRC \n',i);
% material "BeamMat1" is defined below
fprintf(file_inp, '%2.7e \n', r_beam); % geometry data
fprintf(file_inp, '%.6f, %.6f, %.6f \n',n1_elem); % 1-axis of beam element

% Beam Elastic Material
fprintf(file_inp, '*MATERIAL, name=BeamMatl \n');
fprintf(file_inp, '*Elastic \n');
fprintf(file_inp, '%2.5e, %2.5e \n', param_beam);
% param_beam contains two parameter: stiffness and poisson's ratio

% Define Step
fprintf(file_inp, '*Step, name=Step-%d, nlgeom=YES \n',i);
% define steps,  nlgeom=YES option specifies that nonlinear geometry effects
% will be included in the analysis.
fprintf(file_inp, '*Static, STABILIZE=%.4e\n',param_step.STBLfac);
%     fprintf(file_inp, '*Static, STABILIZE\n'); % use automatic stabilization
% specifies a stabilization factor for the static analysis, between [0.1,0.5]
fprintf(file_inp, '%.2e, 1., %.2e, %.2e \n', [param_step.InitStep,param_step.minStep,param_step.maxStep]);
% specifies the initial step size, minimum step size, and maximum step size for the analysis
% "1.' Time period of the step

% ----- Boundary Conditions -----------------------------------
fprintf(file_inp, '*Boundary \n');
% ---- X displacement ----------------------------------------
fprintf(file_inp, 'leftNSET, 1,1,%.3e \n',  disp_BC.Xleft); % apply displacement at left boundary node set
fprintf(file_inp, 'leftNSET, 2,2\n');  % constraint degrees of freedom
fprintf(file_inp, 'leftNSET, 6,6\n');
fprintf(file_inp, 'rightNSET, 1,1,%.3e \n',  disp_BC.Xright);
fprintf(file_inp, 'rightNSET, 2,2\n');
fprintf(file_inp, 'rightNSET, 6,6\n');
% -- Y displacement  -----------------------------------------
fprintf(file_inp, 'lowerNSET, 2,2,%.3e \n',  disp_BC.Ylower);
fprintf(file_inp, 'lowerNSET, 1,1\n');
fprintf(file_inp, 'lowerNSET, 6,6\n');
fprintf(file_inp, 'upperNSET, 2,2,%.3e \n',  disp_BC.Yupper);
fprintf(file_inp, 'upperNSET, 1,1\n');
fprintf(file_inp, 'upperNSET, 6,6\n');
% <node number>, <first dof>, <last dof>, <magnitude of displacement>
% 1,1 option indicates that the node is constrained in one direction only,
% the third field should be blank or equal to the second field.
% Boundary conditions on a node are cumulative
%     fprintf(file_inp, 'FreeRefNode, 2,6\n');
% the node cannot move in degree 2 and 6.
% If the magnitude of the displacement is not specified on the data line, it is assumed to be zero

% Define output
fprintf(file_inp, '*Restart, write, frequency=0 \n');
fprintf(file_inp, '*Output, history, variable=PRESELECT\n');
fprintf(file_inp, '*Output, field, variable=PRESELECT\n');
fprintf(file_inp,'*NODE OUTPUT\n');
fprintf(file_inp,'U, RF,COORD\n');
fprintf(file_inp, '*ELEMENT Output\n');
fprintf(file_inp, 'S,SE,SF,SM,SK,COORD\n');
% S (stress), SE (strain energy), SF (shear force), SM (bending moment), and SK (curvature) will be saved
fprintf(file_inp, '*End Step \n');
% signals the end of the analysis step.
end