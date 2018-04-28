%% global settings
file_dir = fullfile(fileparts(mfilename('fullpath')));

% marker
% markers = {'+' , 'o' , '*' , 'x' , 'square' , 'diamond' , 'v' , '^' , '>' , '<' , 'pentagram' , 'hexagram' , '.', 'none'};
markers = {'+' , 'o' , '*' , 'x' , 'square' , 'diamond' , 'v' , '^' , '>' , '<' , 'pentagram' , 'hexagram' , '.'};
spe_idx = '10';
atom_f = 'HA4';
spe_name = 'OH';
tau = 0.777660157519;
end_t = '0.9';
% end_t = '0.12859156975';
% cycle = 'all';

%% load trajectory data
file_dir2 = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', '..', '..', '..', '..', '..', '..', 'SOHR_DATA');

%% import time
fn_time = fullfile(file_dir2, 'output', 'time_dlsode_M.csv');
delimiter = '';
formatSpec = '%f%[^\n\r]';
%% Open the text file.
fileID = fopen(fn_time,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);
time_vec = dataArray{:, 1};
%% Clear temporary variables
clearvars fn_time delimiter formatSpec fileID dataArray ans;

%% import temperature
fn_temp = fullfile(file_dir2, 'output', 'temperature_dlsode_M.csv');
delimiter = '';
formatSpec = '%f%[^\n\r]';
%% Open the text file.
fileID = fopen(fn_temp,'r');
%% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
%% Close the text file.
fclose(fileID);
%% Allocate imported array to column variable names
temp_vec = dataArray{:, 1};
%% Clear temporary variables
clearvars fn_temp delimiter formatSpec fileID dataArray ans;

%% import time
fn_R = fullfile(file_dir2, 'output', 'reaction_rate_dlsode_M.csv');
delimiter = ',';
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(fn_R,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
%% Close the text file.
fclose(fileID);
%% Create output variable
reaction_R_mat = [dataArray{1:end-1}];
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

%% load pathway probability data

% cycle = 'p1';
% path_idx = [1];
% % delta_t_vec = [1e-5, 5e-4, 1e-3, 2.5e-3, 5e-3, 1e-2, 1.5e-2]; 
% delta_t_vec = [1.2859087486111409e-06, 0.0001285908748611141, 0.001285908748611141, 0.003214771871527852, 0.006429543743055705, 0.009644315614583556];
% y_label_str = '\chi_1';
% title_p = '\color{red} O_2QOOH_1\color{black}\rightarrow(\color{blue}OH\color{black}+)\color{red}OQ^{\prime}OOH_1';

% cycle = 'p2';
% path_idx = [2];
% delta_t_vec = [1.2859087486111409e-06, 0.001285908748611141, 0.01285908748611141, 0.03214771871527852, 0.06429543743055705, 0.09644315614583557];
% y_label_str = '\chi_2';
% title_p = ['\color{red} O_2QOOH_1\color{black}\rightarrow(OH+)\color{red}OQ^{\prime}OOH_1', ... 
%     '\color{black}\rightarrow(\color{blue}OH\color{black}+)\color{red}OQ^{\prime}O_1'];

% cycle = 'p3';
% path_idx = [3];
% delta_t_vec = [1.2859087486111409e-06, 0.001285908748611141, 0.01285908748611141, 0.03214771871527852, 0.06429543743055705, 0.09644315614583557];
% y_label_str = '\chi_3';
% title_p = ['\color{red} O_2QOOH_1\color{black}\rightarrow(OH+)\color{red}OQ^{\prime}OOH_1', ...
%     '\color{black}\rightarrow(OH+)\color{red}OQ^{\prime}O_1' char(10), ...
%     '\color{black}\rightarrow(CH_2O+)\color{red}vinoxy\color{black}', ...
%     '', '[+O_2]\rightarrowCH_2O+CO+\color{blue}OH'];

cycle = 'OHtop100';
path_idx = linspace(1, 100, 100);
% delta_t_vec = [1.2859087486111409e-06, 0.001285908748611141, 0.01285908748611141, 0.1285908748611141, 0.3214771871527852, 0.6429543743055705];
% delta_t_vec = [1.2859087486111409e-06, 0.001285908748611141, 0.01285908748611141, 0.03214771871527852, 0.06429543743055705, 0.09644315614583557];
delta_t_vec = [1.2859087486111409e-06, 0.0001285908748611141, 0.001285908748611141, 0.01285908748611141, 0.03214771871527852, 0.06429543743055705, 0.3214771871527852, 0.6429543743055705];
y_label_str = '$\sum_{i=1}^{100}{\gamma_{i}^{OH}}$';
title_p = ['OH'];

% cycle = 'primary_cycle';
% path_idx = [1, 2, 3];
% delta_t_vec = [0.009644315614583557, 0.09644315614583557];
% % delta_t_vec = [1e-5, 1e-3, 1e-2, 2.5e-2, 5e-2, 1e-1, 2.5e-1];
% y_label_str = '\chi';
% title_p = ['primary cycle'];

% cycle = 'others';

n_path = 100;

% fn_2d_f = fullfile(file_dir, ['Merchant_f_2d_S', spe_idx, '_', atom_f, '_', end_t ,'.csv']);
% fn_2d_f = fullfile(file_dir, ['Merchant_f_2d', '.csv']);
fn_2d_f = fullfile(file_dir, ['Merchant_f_2d_S10_HA4_0.9_100000000_10000_100', '.csv']);
% fn_2d_f = fullfile(file_dir, ['Merchant_f_2d_S61_C_0.9_100000000_10000_100', '.csv']);

delimiter = ',';
formatStr = '%f%f%f';
for i=1:n_path
    formatStr = strcat(formatStr, '%f');
end
formatStr = strcat(formatStr, '%[^\n\r]');
formatSpec = char(formatStr);

%% Open the text file.
fileID = fopen(fn_2d_f,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue', NaN,  'ReturnOnError', false);
%% Close the text file.
fclose(fileID);
f_mat = [dataArray{1:end-1}];

t0 = f_mat(:, 1);
tf = f_mat(:, 2);
% f_value = f_mat(:, end);

tf = tf - t0;

% path index
offset = 2;
% path_idx = linspace(1, n_path, n_path-1+1);
% path_idx = linspace(3, n_path, n_path-3+1);
for i = 1:length(path_idx)
    if i==1
        f_value = f_mat(:, offset + path_idx(i));
    else
        f_value = f_value + f_mat(:, offset + path_idx(i));
    end    
end

% construct 3d surface
xlin = linspace(min(t0), max(t0), 25);
ylin = linspace(min(tf), max(tf), 25);
[X,Y] = meshgrid(xlin, ylin);
f = scatteredInterpolant(t0, tf, f_value);
Z = f(X,Y);

%% update Z
for i = 1:length(X)    
    for j = length(X) - i + 1 : length(X)
        Z(i,j) = nan;
    end
end

%% plot
fig = figure();

%% plot
X_tmp1 = X(1, :);

N = length(delta_t_vec);
colors = lines(N+1);
% colors = colorcube(N);
str_name = cell(N+1,1);
for i=1:N
    str_name{i, 1} = strcat('$\delta$t=', num2str(delta_t_vec(i)*tau,'%1.1e\n'));
end
str_name{N+1, 1} = 'Merchant 3$\times \alpha\beta$';

H = gobjects(N);

for idx=1:N
    % delta
    delta_t = ones(1, length(X_tmp1));
    delta_t = delta_t.* delta_t_vec(idx);
    Z_tmp1 = f(X_tmp1, delta_t);
    % check data
    for i=1:length(X_tmp1)
        if X_tmp1(i) + delta_t(i) > str2double(end_t)
            Z_tmp1(i) = nan;
        end
    end
    H(idx) = plot(X_tmp1 * tau, Z_tmp1, 'LineWidth', 2, 'color', colors(idx, :), ...
        'marker', markers{1, mod(idx-1, length(markers))+ 1}); hold on;
    hold on;
end

idx = idx + 1;
start_idx = 15;
delta_n = 250;
h_beta = plot(time_vec(start_idx:end), 3*alpha(start_idx:end) .* beta(start_idx:end), ...
    'LineWidth', 2, 'color', colors(idx, :));
% h_beta = plot(time_vec(start_idx:end), 3*alpha(start_idx:end) .* beta(start_idx:end), ...
%     'LineWidth', 2, 'color', colors(idx, :), 'HandleVisibility','off');
% hold on;
% scatter(time_vec(start_idx:delta_n:end), 3*alpha(start_idx:delta_n:end) .* beta(start_idx:delta_n:end),...
%         'LineWidth', 2, 'MarkerEdgeColor', colors(idx, :), 'marker', markers{1, mod(idx-1, length(markers))+ 1}, ...
%         'HandleVisibility','off');
% plot(nan, nan, 'LineWidth', 2, 'color', colors(idx, :), 'marker', markers{1, mod(idx-1, length(markers))+ 1});
% hold on;


%% settings
set(gca,'GridLineStyle','--');
xlabel('$t$ (seconds)', 'Interpreter','latex', 'FontSize', 20);
ylabel(y_label_str, 'Interpreter','latex', 'FontSize', 20);
% ylabel(y_label_str, 'FontSize', 20);
% xlim([0, tau*str2double(end_t)]);
xlim([0, max(X_tmp1 * tau)]);
% ylim([0, tau*str2double(end_t)]);

title(title_p, 'Interpreter', 'Latex', 'FontSize', 20);
grid on;

%%  legend
leg_h = legend(str_name, 'Interpreter','latex');
set(leg_h, 'FontSize', 12, 'Box', 'off');
% set(leg_h, 'Location', 'South');

% %% text
% a_x = gca;
% t_x = a_x.XLim(1) + 0.525*(a_x.XLim(2) - a_x.XLim(1));
% t_y = a_x.YLim(1) + 0.278*(a_x.YLim(2) - a_x.YLim(1));
% text(t_x, t_y, [spe_name, '@ $t$' char(10) 'stop path@ $t_f$' char(10) '$t_f > t$'], 'Interpreter','latex', 'FontSize', 20);

%% save to file
figname = strcat('2d_Merchant_f_', end_t, '_S', spe_idx, '_', cycle, '_v4_delta_t_v3.png');
print(fig, fullfile(file_dir, figname), '-r200', '-dpng');


