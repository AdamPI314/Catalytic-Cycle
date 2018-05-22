%% global settings
file_dir = fullfile(fileparts(mfilename('fullpath')));

% marker
% markers = {'+' , 'o' , '*' , 'x' , 'square' , 'diamond' , 'v' , '^' , '>' , '<' , 'pentagram' , 'hexagram' , '.', 'none'};
markers = {'+' , 'o' , '*' , 'x' , 'square' , 'diamond' , 'v' , '^' , '>' , '<' , 'pentagram' , 'hexagram' , '.'};
atom_f = 'HA4';
spe_name = 'npropyl';
tau = 0.777660157519;
end_t = '0.9';

%% global propertities
N = 3;
% delta_t_vec = [1.2859087486111409e-06, 0.001285908748611141, 0.01285908748611141, 0.1285908748611141, 0.3214771871527852, 0.6429543743055705];

% delta_t_value = 0.00001285908748611141;
delta_t_value = 0.0001285908748611141;
% delta_t_value = 0.001285908748611141;
% delta_t_value = 0.009644315614583557;
% delta_t_value = 0.01285908748611141;
% delta_t_value = 0.09644315614583557;
% delta_t_value = 0.1285908748611141;
% delta_t_value = 0.3214771871527852;
% delta_t_value = 0.6429543743055705;

y_label_str = '$\chi$';
time_cell = cell(3, 1);
yvalue_cell = cell(3, 1);

%% load trajectory data
file_dir2 = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', '..', '..', '..', '..', '..', '..', 'SOHR_DATA');

%% import time
fn_time = fullfile(file_dir2, 'output', 'time_dlsode_M.csv');
delimiter3 = '';
formatSpec3 = '%f%[^\n\r]';
%% Open the text file.
fileID3 = fopen(fn_time,'r');
dataArray3 = textscan(fileID3, formatSpec3, 'Delimiter', delimiter3, 'EmptyValue' ,NaN, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID3);
time_vec = dataArray3{:, 1};
%% Clear temporary variables
clearvars fn_time delimiter formatSpec fileID dataArray ans;

%% import temperature
fn_temp = fullfile(file_dir2, 'output', 'temperature_dlsode_M.csv');
delimiter3 = '';
formatSpec3 = '%f%[^\n\r]';
%% Open the text file.
fileID3 = fopen(fn_temp,'r');
%% Read columns of data according to format string.
dataArray3 = textscan(fileID3, formatSpec3, 'Delimiter', delimiter3, 'EmptyValue' ,NaN, 'ReturnOnError', false);
%% Close the text file.
fclose(fileID3);
%% Allocate imported array to column variable names
temp_vec = dataArray3{:, 1};
%% Clear temporary variables
clearvars fn_temp delimiter formatSpec fileID dataArray ans;

%% import time
fn_R = fullfile(file_dir2, 'output', 'reaction_rate_dlsode_M.csv');
delimiter3 = ',';
% For more information, see the TEXTSCAN documentation.
formatSpec3 = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID3 = fopen(fn_R,'r');
dataArray3 = textscan(fileID3, formatSpec3, 'Delimiter', delimiter3, 'EmptyValue' ,NaN, 'ReturnOnError', false);
%% Close the text file.
fclose(fileID3);
%% Create output variable
reaction_R_mat = [dataArray3{1:end-1}];
%% Clear temporary variables
clearvars fn_R delimiter formatSpec fileID dataArray ans;

%% calculate theta first, using equation (13)
r_idx_16 = 1162 + 1;
r_idx_14 = 1080 + 1;
theta = reaction_R_mat(:, r_idx_16) ./ reaction_R_mat(:, r_idx_14);

%% calculate alpha using equation (28)
r_idx_3 = 736 + 1;
r_idx_4 = 738 + 1;
r_idx_20 = 90 + 1;
r_idx_23 = 44 + 1;

alpha = (reaction_R_mat(:, r_idx_3)) ./ (reaction_R_mat(:, r_idx_3) ...
    + reaction_R_mat(:, r_idx_4) ...
    + reaction_R_mat(:, r_idx_20) ...
    + reaction_R_mat(:, r_idx_23));
%% calculate beta using equation (29)
r_idx_12 = 1082 + 1;
r_idx_26 = 914 + 1;
r_idx_27 = 922 + 1;

beta = (theta .*  reaction_R_mat(:, r_idx_14) ) ./ (theta .* reaction_R_mat(:, r_idx_14) ...
    + reaction_R_mat(:, r_idx_12) ...
    + reaction_R_mat(:, r_idx_26) ...
    + reaction_R_mat(:, r_idx_27));

%% read Merchant f value, construct data structure
n_path1 = 100;
fn_2d_f1 = fullfile(file_dir, ['Merchant_f_2d_S60_HA4_0.9_100000000_10000_100', '.csv']);

delimiter1 = ',';
formatStr1 = '%f%f%f';
for i=1:n_path1
    formatStr1 = strcat(formatStr1, '%f');
end
formatStr1 = strcat(formatStr1, '%[^\n\r]');
formatSpec1 = char(formatStr1);

%% Open the text file.
fileID1 = fopen(fn_2d_f1,'r');
dataArray1 = textscan(fileID1, formatSpec1, 'Delimiter', delimiter1, 'EmptyValue', NaN,  'ReturnOnError', false);
%% Close the text file.
fclose(fileID1);
f_mat1 = [dataArray1{1:end-1}];

t01 = f_mat1(:, 1);
tf1 = f_mat1(:, 2);
tf1 = tf1 - t01;

% path index
offset1 = 2;

% nR primary cycle
path_idx1 = linspace(1, 3, 3);
for i = 1:length(path_idx1)
    if i==1
        f_value1 = f_mat1(:, offset1 + path_idx1(i));
    else
        f_value1 = f_value1 + f_mat1(:, offset1 + path_idx1(i));
    end    
end
% construct 3d surface
xlin1 = linspace(min(t01), max(t01), 25);
ylin1 = linspace(min(tf1), max(tf1), 25);
[X1,Y1] = meshgrid(xlin1, ylin1);
f1 = scatteredInterpolant(t01, tf1, f_value1);
Z1 = f1(X1,Y1);
%% update Z
for i = 1:length(X1)    
    for j = length(X1) - i + 1 : length(X1)
        Z1(i,j) = nan;
    end
end
% delta
X_tmp1 = X1(1, :);
delta_t1 = ones(1, length(X_tmp1));
delta_t1 = delta_t1.* delta_t_value;
Z_tmp1 = f1(X_tmp1, delta_t1);
% check data
for i=1:length(X_tmp1)
    if X_tmp1(i) + delta_t1(i) > str2double(end_t)
        Z_tmp1(i) = nan;
    end
end
time_cell{1, 1} = X_tmp1;
yvalue_cell{1, 1} = Z_tmp1;

% nR spur cycle
path_idx2 = linspace(4, n_path1, n_path1-4+1);
for i = 1:length(path_idx2)
    if i==1
        f_value2 = f_mat1(:, offset1 + path_idx2(i));
    else
        f_value2 = f_value2 + f_mat1(:, offset1 + path_idx2(i));
    end    
end
% construct 3d surface
xlin2 = linspace(min(t01), max(t01), 25);
ylin2 = linspace(min(tf1), max(tf1), 25);
[X2,Y2] = meshgrid(xlin2, ylin2);
f2 = scatteredInterpolant(t01, tf1, f_value2);
Z2 = f2(X2,Y2);
%% update Z
for i = 1:length(X2)    
    for j = length(X2) - i + 1 : length(X2)
        Z2(i,j) = nan;
    end
end
% delta
X_tmp2 = X2(1, :);
delta_t2 = ones(1, length(X_tmp2));
delta_t2 = delta_t2.* delta_t_value;
Z_tmp2 = f2(X_tmp2, delta_t2);
% check data
for i=1:length(X_tmp2)
    if X_tmp2(i) + delta_t2(i) > str2double(end_t)
        Z_tmp2(i) = nan;
    end
end
time_cell{2, 1} = X_tmp2;
yvalue_cell{2, 1} = Z_tmp2;

% iR cycle
n_path3 = 100;
fn_2d_f3 = fullfile(file_dir, ['Merchant_f_2d_S61_HA4_0.9_100000000_10000_100', '.csv']);

delimiter3 = ',';
formatStr3 = '%f%f%f';
for i=1:n_path3
    formatStr3 = strcat(formatStr3, '%f');
end
formatStr3 = strcat(formatStr3, '%[^\n\r]');
formatSpec3 = char(formatStr3);

%% Open the text file.
fileID3 = fopen(fn_2d_f3,'r');
dataArray3 = textscan(fileID3, formatSpec3, 'Delimiter', delimiter3, 'EmptyValue', NaN,  'ReturnOnError', false);
%% Close the text file.
fclose(fileID3);
f_mat3 = [dataArray3{1:end-1}];

t03 = f_mat3(:, 1);
tf3 = f_mat3(:, 2);
tf3 = tf3 - t03;

% path index
offset3 = 2;

% nR primary cycle
path_idx3 = linspace(1, n_path3, n_path3);
for i = 1:length(path_idx3)
    if i==1
        f_value3 = f_mat3(:, offset3 + path_idx3(i));
    else
        f_value3 = f_value3 + f_mat3(:, offset3 + path_idx3(i));
    end    
end
% construct 3d surface
xlin3 = linspace(min(t03), max(t03), 25);
ylin3 = linspace(min(tf3), max(tf3), 25);
[X3,Y3] = meshgrid(xlin3, ylin3);
f3 = scatteredInterpolant(t03, tf3, f_value3);
Z3 = f3(X3,Y3);
%% update Z
for i = 1:length(X3)    
    for j = length(X3) - i + 1 : length(X3)
        Z3(i,j) = nan;
    end
end
% delta
X_tmp3 = X3(1, :);
delta_t3 = ones(1, length(X_tmp3));
delta_t3 = delta_t3.* delta_t_value;
Z_tmp3 = f3(X_tmp3, delta_t3);
% check data
for i=1:length(X_tmp3)
    if X_tmp3(i) + delta_t3(i) > str2double(end_t)
        Z_tmp3(i) = nan;
    end
end
time_cell{3, 1} = X_tmp3;
yvalue_cell{3, 1} = Z_tmp3;

%% plot
fig = figure();

%% plot
% colors = lines(N);
% colors = colorcube(N);
colors = [
1.0000         0         0
0         0    1.0000
0.9531    0.6055    0.2578
];

str_name = cell(N,1);
str_name{1, 1} = 'nR primary cycle';
str_name{2, 1} = 'nR spur cycle';
str_name{3, 1} = 'iR cycle';

% stack manually
for idx=2:N
    yvalue_cell{idx,1} = yvalue_cell{idx,1} + yvalue_cell{idx-1,1};
end
% check nan values
for idx=1:N
    validIndices = ~isnan(yvalue_cell{idx,1});
    yvalue_cell{idx,1} = yvalue_cell{idx, 1}(validIndices);
end

validSize = length(yvalue_cell{1,1});
for idx=1:N
    if idx == 1
        polygon_x = [time_cell{idx, 1}(1) * tau time_cell{idx, 1}(1:validSize) * tau time_cell{idx, 1}(validSize) * tau];
        polygon_y = [0 yvalue_cell{idx,1}(1:validSize) 0];
    else
        polygon_x = [time_cell{idx, 1}(1:validSize) * tau fliplr(time_cell{idx, 1}(1:validSize)) * tau];
        polygon_y = [yvalue_cell{idx,1} fliplr(yvalue_cell{idx-1,1})];
    end
    fill(polygon_x, polygon_y, colors(idx, :), ...
        'EdgeColor', colors(idx, :), ...
        'FaceAlpha', 0.9, 'EdgeAlpha', 0.9);
    hold on;
end

%% settings
set(gca,'GridLineStyle','--');
xlabel('$t$ (seconds)', 'Interpreter','latex', 'FontSize', 20);
ylabel(y_label_str, 'Interpreter','latex', 'FontSize', 20);
% ylabel(y_label_str, 'FontSize', 20);
xlim([0, max(time_cell{1,1} * tau)]);

grid on;

%  legend
leg_h = legend(str_name, 'Interpreter','latex');
set(leg_h, 'FontSize', 12, 'Box', 'off');
% set(leg_h, 'Location', 'South');

%% text
a_x = gca;
t_x = a_x.XLim(1) + 0.30*(a_x.XLim(2) - a_x.XLim(1));
t_y = a_x.YLim(1) + 0.40*(a_x.YLim(2) - a_x.YLim(1));
text(t_x, t_y, ['$\delta t=$', num2str(delta_t_value*tau), ' seconds'], 'Interpreter','latex', 'FontSize', 20);

%% save to file
figname = strcat('2d_Merchant_f_', end_t, '_three_cycles_fill_', num2str(delta_t_value*tau), '.png');
print(fig, fullfile(file_dir, figname), '-r200', '-dpng');


