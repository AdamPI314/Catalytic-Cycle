%% global settings
% spe_idx=14;
% spe_name='CO';
spe_idx=17;
spe_name='CH_2O';
tau = '0.777660157519';

%% Current file directory
file_dir = fullfile(fileparts(mfilename('fullpath')));
fig = figure();
n_xaxis = 3;
x_ratio = 1.0;
y_ratio = 1.0/3.0;
%##########################################################################
% Panel 1
%##########################################################################
iax = 1; % Or whichever
subaxis(n_xaxis, 1, iax, 'SpacingVertical',0.05,'SpacingHorizontal',0, ...
    'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0.01,'PaddingBottom',0, ...
    'MarginLeft',.14,'MarginRight',.13,'MarginTop',0.01,'MarginBottom',0.04);
% subplot(3, 1, 1);

folder1 = '0.2seconds';
end_t1 = '0.25718313951098054';
n_path1 = '100';
species_path1 = '';
%% index, delta, every settings
end_idx1_1 = 10;
delta1_1 = 1;
end_idx1_2 = 5;
delta1_2 = 1;
end_idx1_3 = 10;
delta1_3 = 1;

%% const concentration at a time
fn_conc_c1 = fullfile(file_dir, folder1, strcat(species_path1, 'const_conc_', ...
    n_path1, '_',num2str(spe_idx), '_', tau, '_', end_t1, '.csv'));
delimiter = '';
formatSpec = '%f%[^\n\r]';
%% Open the text file.
fileID = fopen(fn_conc_c1,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
%% Close the text file.
fclose(fileID);
%% Allocate imported array to column variable names
const_c1 = dataArray{:, 1};
%% Clear temporary variables
clearvars fn_conc_c1 delimiter formatSpec fileID dataArray ans;

%% sorted pathway probabilities
file_dir = fullfile(fileparts(mfilename('fullpath')));
fn_path1_p1 = fullfile(file_dir, folder1, strcat(species_path1, 'path_prob_top_n_sorted_', ... 
    n_path1, '_',num2str(spe_idx), '_', tau, '_', end_t1, '.csv'));

delimiter = '';
formatSpec = '%f%[^\n\r]';
%% Open the text file.
fileID = fopen(fn_path1_p1,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
%% Close the text file.
fclose(fileID);
path_p_vec1 = dataArray{:, 1};
%% Clear temporary variables
clearvars fn_path1_p1 delimiter formatSpec fileID dataArray ans;

%% cumulative pathway probabilities
cumu_path_p_vec1 = path_p_vec1;
for i = 2:length(cumu_path_p_vec1)
    cumu_path_p_vec1(i) = cumu_path_p_vec1(i) + cumu_path_p_vec1(i-1);
end

%% linear array
path_idx1 = linspace(1, length(path_p_vec1), length(path_p_vec1));
path_idx1 = path_idx1';

%% plot
% fig = figure();
% https://www.mathworks.com/help/matlab/graphics_transition/why-are-plot-lines-different-colors.html
% https://www.mathworks.com/help/matlab/creating_plots/customize-graph-with-two-y-axes.html
co = [    0    0.4470    0.7410 % 1th plot
%     0.8500    0.3250    0.0980 % 2nd plot
    0.9290    0.6940    0.1250 % 3rd plot
%     0.4940    0.1840    0.5560 % 4th plot
%     0.4660    0.6740    0.1880 % 5th plot
%     0.3010    0.7450    0.9330 % 6th plot
%     0.6350    0.0780    0.1840 % 7th plot
%     0         0    1.0000 % 8th plot, blue
    1   0   0 % for placeholder
    1   0   0]; % 9th plot, red
set(fig,'defaultAxesColorOrder',co)

% plot pathway probabilities
% individual
p1 = plot([path_idx1(1), path_idx1(end_idx1_1)], [const_c1, const_c1], ...
    'color',co(1, :), 'LineWidth', 2); hold on;
% cumulative
p2 = plot(path_idx1(1:delta1_1:end_idx1_1), cumu_path_p_vec1(1:delta1_1:end_idx1_1) , ...
    'color',co(2, :), 'linestyle', '-', 'LineWidth', 2); hold on;
% placeholder for legend
p3 = plot(NaN, NaN, 'LineStyle', ':', 'LineWidth', 2, 'color', co(3, :));
%% conc
set(gca,'GridLineStyle','--');
% xlabel('#Path', 'FontSize', 20);
ylabel('Concentration', 'FontSize', 20);
ylim([10^(1.0*log10(cumu_path_p_vec1(1))), 10^(0.9995*(log10(const_c1)))]);
% ylim([0.0, 10^(0.9995*(log10(const_c1)))]);

a_x1 = gca;
t_x1 = a_x1.XLim(1) + 0.250*(a_x1.XLim(2) - a_x1.XLim(1));
t_y1 = a_x1.YLim(1) + 0.850*(a_x1.YLim(2) - a_x1.YLim(1));
text(t_x1, t_y1, [spe_name, '@', 't=stage1A', ', primitive pathways']);

yyaxis right;
error_data1 = (const_c1 - cumu_path_p_vec1)./const_c1.*100;
semilogy(path_idx1(1:delta1_1:end_idx1_1), error_data1(1:delta1_1:end_idx1_1), ... 
    'LineStyle', ':', 'LineWidth', 2, 'color', co(3, :)); hold on;
ylabel('%Relative Error', 'FontSize', 20);
% set(gca, 'ytick', []);

%% global settings
grid on;
xlim([1, end_idx1_1]);
leg_h = legend([p1; p2; p3],'EXACT','SOHR', 'Percentage Error');
set(leg_h, 'FontSize', 14, 'Box', 'off');
% [left, bottom, weight, height]
set(leg_h, 'Position', [0.275*x_ratio, (0.19+n_xaxis-iax)*y_ratio, 0.5*x_ratio, 0.1*y_ratio]);

%% Zoom in figure1
% create a new pair of axes inside current figure
z1_position = [.23*x_ratio (.45+n_xaxis-iax)*y_ratio .25*x_ratio .25*y_ratio];
axes('position',z1_position);
box on; % put box around new pair of axes
plot([path_idx1(1), path_idx1(end_idx1_2)], [const_c1, const_c1], ...
    'color',co(1, :), 'LineWidth', 2); hold on;
% cumulative
semilogy(path_idx1(1:delta1_2:end_idx1_2), cumu_path_p_vec1(1:delta1_2:end_idx1_2) , ...
    'color',co(2, :), 'linestyle', '-', 'LineWidth', 2); hold on;
axis tight;
grid on;
ylim([10^(1.0*log10(cumu_path_p_vec1(1))), 10^(0.999*(log10(const_c1)))]);

%% Error
yyaxis right;
semilogy(path_idx1(1:delta1_2:end_idx1_2), error_data1(1:delta1_2:end_idx1_2), ... 
    'LineStyle', ':', 'LineWidth', 2, 'color', 'r'); hold on;

%% Zoom in figure2
co = [    0    0.4470    0.7410 % 1th plot
%     0.8500    0.3250    0.0980 % 2nd plot
%       0.9290    0.6940    0.1250 % 3rd plot
    0.4940    0.1840    0.5560 % 4th plot
%     0.4660    0.6740    0.1880 % 5th plot
%       0.3010    0.7450    0.9330 % 6th plot
%     0.6350    0.0780    0.1840 % 7th plot
%     0         0    1.0000 % 8th plot, blue
    1   0   0 % for placeholder
    1   0   0]; % 9th plot, red
set(fig,'defaultAxesColorOrder',co);
z2_position = [.55*x_ratio (.45+n_xaxis-iax)*y_ratio .25*x_ratio .25*y_ratio];
% create a new pair of axes inside current figure
axes('position',z2_position);

% plot a placeholder
plot(NaN, NaN);

set(gca, 'ytick', []);
yyaxis right
box on; % put box around new pair of axes
% cumulative
zf2_1 = semilogy(path_idx1(1:delta1_3:end_idx1_3), path_p_vec1(1:delta1_3:end_idx1_3) , ...
    'color',co(2, :), 'linestyle', '-.', 'LineWidth', 2); hold on;
tick_pos = linspace(ceil(log10(path_p_vec1(end_idx1_3))),floor(log10(path_p_vec1(1))), ...
    1+floor(log10(path_p_vec1(1))) - ceil(log10(path_p_vec1(end_idx1_3))));
tick_pos = arrayfun(@(x) 10^x, tick_pos);
set(gca, 'ytick', tick_pos);
axis tight;
grid on;
% legend
leg_z2 = legend(zf2_1,['Pathway' newline 'Probability']);
set(leg_z2, 'FontSize', 10, 'Box', 'off');
% [left, bottom, weight, height]
set(leg_z2, 'Position', [z2_position(1) + z2_position(3)*0.375, ...
    z2_position(2) + z2_position(4)*0.5, 0.1*x_ratio, 0.1*y_ratio]);
% alpha(0.5);

%##########################################################################
% Panel 2
%##########################################################################
iax = 2; % Or whichever
subaxis(n_xaxis, 1, iax, 'SpacingVertical',0.04,'SpacingHorizontal',0, ...
    'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
    'MarginLeft',.14,'MarginRight',.13,'MarginTop',0,'MarginBottom',0.04);
% subplot(3, 1, 2);

folder2 = '0.5tau';
end_t2 = '0.5';
n_path2 = '500';
species_path2 = '';
%% index, delta, every settings
end_idx2_1 = 500;
delta2_1 = 1;
end_idx2_2 = 50;
delta2_2 = 1;
end_idx2_3 = 500;
delta2_3 = 1;

%% const concentration at a time
fn_conc_c1 = fullfile(file_dir, folder2, strcat(species_path2, 'const_conc_', ...
    n_path2, '_',num2str(spe_idx), '_', tau, '_', end_t2, '.csv'));
delimiter = '';
formatSpec = '%f%[^\n\r]';
%% Open the text file.
fileID = fopen(fn_conc_c1,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
%% Close the text file.
fclose(fileID);
%% Allocate imported array to column variable names
const_c1 = dataArray{:, 1};
%% Clear temporary variables
clearvars fn_conc_c1 delimiter formatSpec fileID dataArray ans;

%% sorted pathway probabilities
file_dir = fullfile(fileparts(mfilename('fullpath')));
fn_path2_p1 = fullfile(file_dir, folder2, strcat(species_path2, 'path_prob_top_n_sorted_', ... 
    n_path2, '_',num2str(spe_idx), '_', tau, '_', end_t2, '.csv'));

delimiter = '';
formatSpec = '%f%[^\n\r]';
%% Open the text file.
fileID = fopen(fn_path2_p1,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
%% Close the text file.
fclose(fileID);
path_p_vec2 = dataArray{:, 1};
%% Clear temporary variables
clearvars fn_path2_p1 delimiter formatSpec fileID dataArray ans;

%% cumulative pathway probabilities
cumu_path_p_vec2 = path_p_vec2;
for i = 2:length(cumu_path_p_vec2)
    cumu_path_p_vec2(i) = cumu_path_p_vec2(i) + cumu_path_p_vec2(i-1);
end

%% linear array
path_idx2 = linspace(1, length(path_p_vec2), length(path_p_vec2));
path_idx2 = path_idx2';

%% plot
% fig = figure();
% https://www.mathworks.com/help/matlab/graphics_transition/why-are-plot-lines-different-colors.html
% https://www.mathworks.com/help/matlab/creating_plots/customize-graph-with-two-y-axes.html
co = [    0    0.4470    0.7410 % 1th plot
%     0.8500    0.3250    0.0980 % 2nd plot
    0.9290    0.6940    0.1250 % 3rd plot
%     0.4940    0.1840    0.5560 % 4th plot
%     0.4660    0.6740    0.1880 % 5th plot
%     0.3010    0.7450    0.9330 % 6th plot
%     0.6350    0.0780    0.1840 % 7th plot
%     0         0    1.0000 % 8th plot, blue
    1   0   0 % for placeholder
    1   0   0]; % 9th plot, red
set(fig,'defaultAxesColorOrder',co)

% plot pathway probabilities
% individual
p1 = plot([path_idx2(1), path_idx2(end_idx2_1)], [const_c1, const_c1], ...
    'color',co(1, :), 'LineWidth', 2); hold on;
% cumulative
p2 = plot(path_idx2(1:delta2_1:end_idx2_1), cumu_path_p_vec2(1:delta2_1:end_idx2_1) , ...
    'color',co(2, :), 'linestyle', '-', 'LineWidth', 2); hold on;
% placeholder for legend
p3 = plot(NaN, NaN, 'LineStyle', ':', 'LineWidth', 2, 'color', co(3, :));
%% conc
set(gca,'GridLineStyle','--');
% xlabel('#Path', 'FontSize', 20);
ylabel('Concentration', 'FontSize', 20);
ylim([10^(1.0*log10(cumu_path_p_vec2(1))), 10^(0.9995*(log10(const_c1)))]);
% ylim([0.0, 10^(0.9995*(log10(const_c1)))]);

a_x2 = gca;
t_x2 = a_x2.XLim(1) + 0.250*(a_x2.XLim(2) - a_x2.XLim(1));
t_y2 = a_x2.YLim(1) + 0.850*(a_x2.YLim(2) - a_x2.YLim(1));
text(t_x2, t_y2, [spe_name, '@', 't=0.5\tau', ', primitive pathways']);

yyaxis right;
error_data1 = (const_c1 - cumu_path_p_vec2)./const_c1.*100;
semilogy(path_idx2(1:delta2_1:end_idx2_1), error_data1(1:delta2_1:end_idx2_1), ... 
    'LineStyle', ':', 'LineWidth', 2, 'color', co(3, :)); hold on;
ylabel('%Relative Error', 'FontSize', 20);
% set(gca, 'ytick', []);

%% global settings
grid on;
xlim([1, end_idx2_1]);
% leg_h = legend([p1; p2; p3],'EXACT','SOHR', 'Percentage Error');
% set(leg_h, 'FontSize', 14, 'Box', 'off');
% % [left, bottom, weight, height]
% set(leg_h, 'Position', [0.275*x_ratio, (0.19+n_xaxis-iax)*y_ratio, 0.5*x_ratio, 0.1*y_ratio]);

%% Zoom in figure1
% create a new pair of axes inside current figure
z1_position = [.23*x_ratio (.45+n_xaxis-iax)*y_ratio .25*x_ratio .25*y_ratio];
axes('position',z1_position);
box on; % put box around new pair of axes
plot([path_idx2(1), path_idx2(end_idx2_2)], [const_c1, const_c1], ...
    'color',co(1, :), 'LineWidth', 2); hold on;
% cumulative
semilogy(path_idx2(1:delta2_2:end_idx2_2), cumu_path_p_vec2(1:delta2_2:end_idx2_2) , ...
    'color',co(2, :), 'linestyle', '-', 'LineWidth', 2); hold on;
axis tight;
grid on;
ylim([10^(1.0*log10(cumu_path_p_vec2(1))), 10^(0.999*(log10(const_c1)))]);

%% Error
yyaxis right
semilogy(path_idx2(1:delta2_2:end_idx2_2), error_data1(1:delta2_2:end_idx2_2), ... 
    'LineStyle', ':', 'LineWidth', 2, 'color', 'r'); hold on;


%% Zoom in figure2
co = [    0    0.4470    0.7410 % 1th plot
%     0.8500    0.3250    0.0980 % 2nd plot
%       0.9290    0.6940    0.1250 % 3rd plot
    0.4940    0.1840    0.5560 % 4th plot
%     0.4660    0.6740    0.1880 % 5th plot
%       0.3010    0.7450    0.9330 % 6th plot
%     0.6350    0.0780    0.1840 % 7th plot
%     0         0    1.0000 % 8th plot, blue
    1   0   0 % for placeholder
    1   0   0]; % 9th plot, red
set(fig,'defaultAxesColorOrder',co);
z2_position = [.55*x_ratio (.45+n_xaxis-iax)*y_ratio .25*x_ratio .25*y_ratio];
% create a new pair of axes inside current figure
axes('position',z2_position);

% plot a placeholder
plot(NaN, NaN);

set(gca, 'ytick', []);
yyaxis right
box on; % put box around new pair of axes
% cumulative
zf2_1 = semilogy(path_idx2(1:delta2_3:end_idx2_3), path_p_vec2(1:delta2_3:end_idx2_3) , ...
    'color',co(2, :), 'linestyle', '-.', 'LineWidth', 2); hold on;
tick_pos = linspace(ceil(log10(path_p_vec2(end_idx2_3))),floor(log10(path_p_vec2(1))), ...
    1+floor(log10(path_p_vec2(1))) - ceil(log10(path_p_vec2(end_idx2_3))));
tick_pos = arrayfun(@(x) 10^x, tick_pos);
set(gca, 'ytick', tick_pos);
axis tight;
grid on;
% legend
leg_z2 = legend(zf2_1,['Pathway' newline 'Probability']);
set(leg_z2, 'FontSize', 10, 'Box', 'off');
% [left, bottom, weight, height]
set(leg_z2, 'Position', [z2_position(1) + z2_position(3)*0.375, ...
    z2_position(2) + z2_position(4)*0.5, 0.1*x_ratio, 0.1*y_ratio]);
% alpha(0.5);

%##########################################################################
% Panel 3
%##########################################################################
iax = 3; % Or whichever
subaxis(n_xaxis, 1, iax, 'SpacingVertical',0.06,'SpacingHorizontal',0, ...
    'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0.0, ...
    'MarginLeft',.14,'MarginRight',.13,'MarginTop',0.0,'MarginBottom',0.06);
% subplot(3, 1, 3);

folder3 = '0.5tau2';
end_t3 = '0.5';
n_path3 = '50';
species_path3 = 'species_';
%% index, delta, every settings
end_idx3_1 = 50;
delta3_1 = 1;
end_idx3_2 = 10;
delta3_2 = 1;
end_idx3_3 = 50;
delta3_3 = 1;

%% const concentration at a time
fn_conc_c1 = fullfile(file_dir, folder3, strcat(species_path3, 'const_conc_', ...
    n_path3, '_',num2str(spe_idx), '_', tau, '_', end_t3, '.csv'));
delimiter = '';
formatSpec = '%f%[^\n\r]';
%% Open the text file.
fileID = fopen(fn_conc_c1,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
%% Close the text file.
fclose(fileID);
%% Allocate imported array to column variable names
const_c1 = dataArray{:, 1};
%% Clear temporary variables
clearvars fn_conc_c1 delimiter formatSpec fileID dataArray ans;

%% sorted pathway probabilities
file_dir = fullfile(fileparts(mfilename('fullpath')));
fn_path3_p1 = fullfile(file_dir, folder3, strcat(species_path3, 'path_prob_top_n_sorted_', ... 
    n_path3, '_',num2str(spe_idx), '_', tau, '_', end_t3, '.csv'));

delimiter = '';
formatSpec = '%f%[^\n\r]';
%% Open the text file.
fileID = fopen(fn_path3_p1,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
%% Close the text file.
fclose(fileID);
path_p_vec3 = dataArray{:, 1};
%% Clear temporary variables
clearvars fn_path3_p1 delimiter formatSpec fileID dataArray ans;

%% cumulative pathway probabilities
cumu_path_p_vec3 = path_p_vec3;
for i = 2:length(cumu_path_p_vec3)
    cumu_path_p_vec3(i) = cumu_path_p_vec3(i) + cumu_path_p_vec3(i-1);
end

%% linear array
path_idx3 = linspace(1, length(path_p_vec3), length(path_p_vec3));
path_idx3 = path_idx3';

%% plot
% fig = figure();
% https://www.mathworks.com/help/matlab/graphics_transition/why-are-plot-lines-different-colors.html
% https://www.mathworks.com/help/matlab/creating_plots/customize-graph-with-two-y-axes.html
co = [    0    0.4470    0.7410 % 1th plot
%     0.8500    0.3250    0.0980 % 2nd plot
    0.9290    0.6940    0.1250 % 3rd plot
%     0.4940    0.1840    0.5560 % 4th plot
%     0.4660    0.6740    0.1880 % 5th plot
%     0.3010    0.7450    0.9330 % 6th plot
%     0.6350    0.0780    0.1840 % 7th plot
%     0         0    1.0000 % 8th plot, blue
    1   0   0 % for placeholder
    1   0   0]; % 9th plot, red
set(fig,'defaultAxesColorOrder',co)

% plot pathway probabilities
% individual
p1 = plot([path_idx3(1), path_idx3(end_idx3_1)], [const_c1, const_c1], ...
    'color',co(1, :), 'LineWidth', 2); hold on;
% cumulative
p2 = plot(path_idx3(1:delta3_1:end_idx3_1), cumu_path_p_vec3(1:delta3_1:end_idx3_1) , ...
    'color',co(2, :), 'linestyle', '-', 'LineWidth', 2); hold on;
% placeholder for legend
p3 = plot(NaN, NaN, 'LineStyle', ':', 'LineWidth', 2, 'color', co(3, :));
%% conc
set(gca,'GridLineStyle','--');
xlabel('Path#', 'FontSize', 20);
ylabel('Concentration', 'FontSize', 20);
ylim([10^(1.0*log10(cumu_path_p_vec3(1))), 10^(0.9995*(log10(const_c1)))]);
% ylim([0.0, 10^(0.9995*(log10(const_c1)))]);

%% subfigure text
a_x3 = gca;
t_x3 = a_x3.XLim(1) + 0.250*(a_x3.XLim(2) - a_x3.XLim(1));
t_y3 = a_x3.YLim(1) + 0.850*(a_x3.YLim(2) - a_x3.YLim(1));
text(t_x3, t_y3, [spe_name, '@', 't=0.5\tau', ', merged pathways']);

yyaxis right;
error_data1 = (const_c1 - cumu_path_p_vec3)./const_c1.*100;
semilogy(path_idx3(1:delta3_1:end_idx3_1), error_data1(1:delta3_1:end_idx3_1), ... 
    'LineStyle', ':', 'LineWidth', 2, 'color', co(3, :)); hold on;
ylabel('%Relative Error', 'FontSize', 20);
% set(gca, 'ytick', []);

%% global settings
grid on;
xlim([1, end_idx3_1]);
% leg_h = legend([p1; p2; p3],'EXACT','SOHR', 'Percentage Error');
% set(leg_h, 'FontSize', 14, 'Box', 'off');
% % [left, bottom, weight, height]
% set(leg_h, 'Position', [0.275*x_ratio, (0.19+n_xaxis-iax)*y_ratio, 0.5*x_ratio, 0.1*y_ratio]);

%% Zoom in figure1
% create a new pair of axes inside current figure
z1_position = [.23*x_ratio (.45+n_xaxis-iax)*y_ratio .25*x_ratio .25*y_ratio];
axes('position',z1_position);
box on; % put box around new pair of axes
plot([path_idx3(1), path_idx3(end_idx3_2)], [const_c1, const_c1], ...
    'color',co(1, :), 'LineWidth', 2); hold on;
% cumulative
semilogy(path_idx3(1:delta3_2:end_idx3_2), cumu_path_p_vec3(1:delta3_2:end_idx3_2) , ...
    'color',co(2, :), 'linestyle', '-', 'LineWidth', 2); hold on;
axis tight;
grid on;
ylim([10^(1.0*log10(cumu_path_p_vec3(1))), 10^(0.999*(log10(const_c1)))]);

%% Error
yyaxis right;
semilogy(path_idx3(1:delta3_2:end_idx3_2), error_data1(1:delta3_2:end_idx3_2), ... 
    'LineStyle', ':', 'LineWidth', 2, 'color', 'r'); hold on;

%% Zoom in figure2
co = [    0    0.4470    0.7410 % 1th plot
%     0.8500    0.3250    0.0980 % 2nd plot
%       0.9290    0.6940    0.1250 % 3rd plot
    0.4940    0.1840    0.5560 % 4th plot
%     0.4660    0.6740    0.1880 % 5th plot
%       0.3010    0.7450    0.9330 % 6th plot
%     0.6350    0.0780    0.1840 % 7th plot
%     0         0    1.0000 % 8th plot, blue
    1   0   0 % for placeholder
    1   0   0]; % 9th plot, red
set(fig,'defaultAxesColorOrder',co);
z2_position = [.55*x_ratio (.45+n_xaxis-iax)*y_ratio .25*x_ratio .25*y_ratio];
% create a new pair of axes inside current figure
axes('position',z2_position);

% plot a placeholder
plot(NaN, NaN);

set(gca, 'ytick', []);
yyaxis right;
box on; % put box around new pair of axes
% cumulative
zf2_1 = semilogy(path_idx3(1:delta3_3:end_idx3_3), path_p_vec3(1:delta3_3:end_idx3_3) , ...
    'color',co(2, :), 'linestyle', '-.', 'LineWidth', 2); hold on;
tick_pos = linspace(ceil(log10(path_p_vec3(end_idx3_3))),floor(log10(path_p_vec3(1))), ...
    1+floor(log10(path_p_vec3(1))) - ceil(log10(path_p_vec3(end_idx3_3))));
tick_pos = arrayfun(@(x) 10^x, tick_pos);
set(gca, 'ytick', tick_pos);
axis tight;
grid on;
% legend
leg_z2 = legend(zf2_1,['Pathway' newline 'Probability']);
set(leg_z2, 'FontSize', 10, 'Box', 'off');
% [left, bottom, weight, height]
set(leg_z2, 'Position', [z2_position(1) + z2_position(3)*0.375, ...
    z2_position(2) + z2_position(4)*0.5, 0.1*x_ratio, 0.1*y_ratio]);
% alpha(0.5);

%##########################################################################
% Figure size
%##########################################################################
x0=10;
y0=10;
width=400;
height=900;
set(gcf,'units','points','position',[x0,y0,width,height]);

%% save to file
figname = strcat(species_path1, 'pathway_prob_concentration_3in1_',num2str(spe_idx), '_',end_t1, '.png');
print(fig, fullfile(file_dir, figname), '-r200', '-dpng');