clear
close all

% --------------------------------------------------------------------------
% Postprocess numerical simulation result of the discrete fiber network
% calculate stress and strain, and compare with homogenized model
% -------------------------------------------------------------------------

%% Preliminary parameter
width = 49; % width of the fiber network
thickness = 2;
volume = width * width * thickness; % volume of the fiber network
num_interp = 20; % number of interval for interpolation
scale_mpa = 1; % scaling factor for stress

%% Read simulation data from xlsx file (discrete fiber network)
% ---- Import initial coordinate of node set -----------------------------
dir_coord = uigetdir(' ','Select initial coordinate data folder');
fulldir_coord = fullfile(dir_coord,'*.mat');
list_coord = dir(fulldir_coord); % get list of data path

% ------ Import equi-displacement U at right nodeset --------------------------------------
dir_equiU = uigetdir(' ','Select equi U displacement data folder');
fulldir_equiU = fullfile(dir_equiU,'*.xlsx');
list_pathequiU = dir(fulldir_equiU); % get list of data path

% ------ Import equi-displacement U1 --------------------------------------
dir_U1 = uigetdir(' ','Select U1 displacement data folder');
fulldir_U1 = fullfile(dir_U1,'*.xlsx');
list_pathU1 = dir(fulldir_U1); % get list of data path

% ------ Import equi-displacement U2 --------------------------------------
dir_U2 =  uigetdir(' ','Select U2 data folder');
fulldir_U2 = fullfile(dir_U2,'*.xlsx');
list_pathU2 = dir(fulldir_U2); % get list of data path

% ------ Import reaction force in x direction: RF1 --------------------------------------
dir_RF1 = uigetdir(' ','Select RF1 data folder');
fulldir_RF1 = fullfile(dir_RF1,'*.xlsx');
list_pathRF1 = dir(fulldir_RF1); % get list of data path

% ------ Import reaction force in y direction: RF2 --------------------------------------
dir_RF2 = uigetdir(' ','Select RF1 data folder');
fulldir_RF2 = fullfile(dir_RF2,'*.xlsx');
list_pathRF2 = dir(fulldir_RF2); % get list of data path

%% Postprocessing of simulation data
sigma_xx = [];
sigma_yy = [];
sigma_xy = [];
sigma_yx = [];

for i = 1:size(list_pathRF1,1)
     fprintf('Postprocess started for discrete fiber network %.f\n',i)
    % ---- calculate strain ------------------------------------------
    name_temp = list_pathequiU(i).name;
    path_temp = fullfile(dir_equiU, name_temp); % construct complete file path
    data_U  = readmatrix(path_temp);
    temp_U = data_U(:,2);
    epsilon = temp_U .* (2 / width * 100);
    epsilon_interp = (linspace(0, max(epsilon), num_interp));
    
    % ----- import U1 and U2 -------------------------------------------
    U1 = fieldImport(list_pathU1(i),dir_U1);
    U2 = fieldImport(list_pathU2(i),dir_U2);
    
    % --- calculate deformed coordinate ----------------------------------
    [coord_x,coord_y] = getDeformed(list_coord(i),dir_coord,U1,U2);
    
    % ---- import RF1 and RF2 --------------------------------------------
    RF1 = fieldImport(list_pathRF1(i),dir_RF1);
    RF2 = fieldImport(list_pathRF2(i),dir_RF2);
    
    % ---- calculate stress ----------------------------------------------
    sigma_xx = getStress(coord_x,RF1,volume,sigma_xx,epsilon,epsilon_interp);
    sigma_yy = getStress(coord_y,RF2,volume,sigma_yy,epsilon,epsilon_interp);
    sigma_xy = getStress(coord_x,RF2,volume,sigma_xy,epsilon,epsilon_interp);
    sigma_yx = getStress(coord_y,RF1,volume,sigma_yx,epsilon,epsilon_interp);
    
    fprintf('Postprocess finished for discrete fiber network %.f\n',i)
end

%% ----- calculate mean and standard deviation ----------------------
[mean_xx, pos_xx, neg_xx,error_xx] = dataPost(epsilon_interp, sigma_xx / scale_mpa,...
    'Sxx',[0.7,0.7,0.7],[0.5,0.5,0.5],'Sxx (MPa)');
[mean_yy, pos_yy, neg_yy,error_yy] = dataPost(epsilon_interp, sigma_yy / scale_mpa,...
    'Syy',[0.7,0.7,0.7],[0.5,0.5,0.5],'Syy (MPa)');
[mean_xy, pos_xy, neg_xy,error_xy] = dataPost(epsilon_interp, sigma_xy / scale_mpa,...
    'Sxy',[1,0.8,0.8],[1,0,0],'Sxy (MPa)');
[mean_yx, pos_yx, neg_yx,error_yx] = dataPost(epsilon_interp, sigma_yx / scale_mpa,...
    'Syx',[1,0.8,0.8],[1,0,0],'Syx (MPa)');

%% ---- compare maximum ------------------
compareMax(mean_xx,mean_yy,error_xx,error_yy,'Sxx','Syy','MaxNormal',[0.5,0.5,0.5],'Normal stress (MPa)')
compareMax(mean_xy,mean_yx,error_xy,error_yx,'Sxy','Syx','Maxshear',[0.5,0,0],'Shear stress (MPa)')

%% ---- Compare with homogenized fiber model (HGO) result ---------------------------------
[file,path] = uigetfile(' ','Select the numerical result of HGO');
dataHGO = readmatrix(fullfile(path,file));
Sxx_hgo = dataHGO(:,2); % homogenized stress
Syy_hgo = dataHGO(:,4);
Sxy_hgo = dataHGO(:,6);
strain_hgo = dataHGO(:,9);

% --- PLot stress-strain curve --------------------------------------------
figure()
plot(strain_hgo,Sxx_hgo,'LineWidth',2,'Color',[0.5,0.5,0.5])
xlabel('Strain(%)')
ylabel('Sxx(MPa)')
ylim([0 8e5])

figure()
plot(strain_hgo,Syy_hgo,'LineWidth',2,'Color',[0.5,0.5,0.5])
xlabel('Strain(%)')
ylabel('Syy(MPa)')
ylim([0 8e5])

figure()
plot(strain_hgo,Sxy_hgo,'LineWidth',2,'Color',[1,0,0])
xlabel('Strain(%)')
ylabel('Sxy(MPa)')
ylim([0 8e5])

% ---- compare maximum -------------------------------------------
max_S = [max(Sxx_hgo) max(Syy_hgo)];
x = categorical({'Sxx','Syy'});
x_ordered = reordercats(x,{'Sxx','Syy'});

f3 = figure('Name','MaxS_HGO','NumberTitle','off');
figure(f3)
bar_plot = bar(x_ordered,max_S);
bar_plot.FaceColor = 'flat';
bar_plot.CData(1,:) = [0.7,0.7,0.7];
bar_plot.CData(2,:) = [0.7,0.7,0.7];

%% Function
function myData = dataImport(list_path,dir_data)
name_temp = list_path.name;
path_temp = fullfile(dir_data, name_temp); % construct complete file path
myData  = readmatrix(path_temp);
end

function U = fieldImport(path_U,dir_U)
data_U = dataImport(path_U,dir_U);
U = [];
for k = 1:2:size(data_U,2) - 1
    temp_U = data_U(:,k+1)'; % a row represent a node
    U = [U;temp_U];
end
end

function [coord_x,coord_y] = getDeformed(path_coord,dir_coord,U1,U2)
% --- import initial coordinate ----------------------------------------
name_temp = path_coord.name;
path_temp = fullfile(dir_coord, name_temp); % construct complete file path
data_coord0  = load(path_temp);
coord0_struct = data_coord0.coord_BC;
coord0_idxy = [coord0_struct.left;coord0_struct.right;coord0_struct.lower;coord0_struct.upper];
[coord0_sorted,idx_sorted] = sortrows(coord0_idxy,1);
coord_x = coord0_sorted(:,2) + U1;
coord_y = coord0_sorted(:,3) + U2;

% % -------- plot to check the deformed coordinates ---------------------------------------------
% figure()
% plot(coord0_sorted(:,2),coord0_sorted(:,3),'b.')
% hold on
% plot(coord_x(:,end),coord_y(:,end),'rx')
% % ---------------------------------------------------------------------
end

function list_sigma = getStress(coord,RF,volume,list_sigma,x_initial,interp_x)
% ---- Stress of discrete fiber network calculated as 
% the volume averaging of microscale stress across the model --------------
temp_sigma = sum(coord .* RF,1) / volume;
interp_sigma = interp1(x_initial, temp_sigma, interp_x); % interpolate the stress
list_sigma = [list_sigma;interp_sigma];
end



function [mean_Y, dev_pos, dev_neg,error_max] = dataPost(X_data, Y_data,fig_name,shade_color,line_color,ylabel_name)
% --- calculate the mean stress and standard deviation -----------
mean_Y = mean(Y_data,1);
std_Y = std(Y_data,0,1);
error_max = std_Y(end);

dev_pos = mean_Y + std_Y; % standard deviation range
dev_neg = mean_Y - std_Y;

% --- plot the mean stress-strain curve, and standard deviation is plotted as shaded region -----------
f1 = figure('Name',fig_name,'NumberTitle','off');
figure(f1)
axes('NextPlot', 'add'); 
patch([X_data fliplr(X_data)], [dev_pos  fliplr(dev_neg)],shade_color,'EdgeColor','none')
hold on
plot(X_data, mean_Y,'Color',line_color);
xlabel('Strain(%)')
ylabel(ylabel_name)
end

function compareMax(mean_S11,mean_S22, error_11,error_22,xlabel_1,xlabel_2,fig_name,bar_color,ylabel_name)
max_S = [mean_S11(end) mean_S22(end)];
x = categorical({xlabel_1,xlabel_2});
x_ordered = reordercats(x,{xlabel_1,xlabel_2});
error_neg = [error_11 error_22];
error_pos = error_neg;

f3 = figure('Name',fig_name,'NumberTitle','off');
figure(f3)
bar_plot = bar(x_ordered,max_S);
bar_plot.FaceColor = 'flat';
bar_plot.CData(1,:) = bar_color;
bar_plot.CData(2,:) = bar_color;
% bar_plot(2).Cdata = [0.7,0.7,0.7];

hold on
er = errorbar(x,max_S,error_neg,error_pos);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
% ylim([0 100])
ylabel(ylabel_name)
end
