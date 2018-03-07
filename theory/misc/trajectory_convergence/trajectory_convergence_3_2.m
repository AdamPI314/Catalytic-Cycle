%% global settings
spe_idx=59;
spe_name='C_3H_6';
folder1 = '0.5tau';
tau = '0.777660157519';
end_t = '0.5';

%% Current file directory
file_dir = fullfile(fileparts(mfilename('fullpath')));

%% const concentration at a time
fn_conc_c1 = fullfile(file_dir, folder1, strcat('const_conc_500_',num2str(spe_idx),'_',tau,'_',end_t,'.csv'));
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
fn_path_p1 = fullfile(file_dir, folder1, strcat('path_prob_top_n_sorted_500_',num2str(spe_idx),'_',tau,'_',end_t,'.csv'));

delimiter = '';
formatSpec = '%f%[^\n\r]';
%% Open the text file.
fileID = fopen(fn_path_p1,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
%% Close the text file.
fclose(fileID);
path_p_vec1 = dataArray{:, 1};
%% Clear temporary variables
clearvars fn_path_p1 delimiter formatSpec fileID dataArray ans;

%% cumulative pathway probabilities
cumu_path_p_vec1 = path_p_vec1;
for i = 2:length(cumu_path_p_vec1)
    cumu_path_p_vec1(i) = cumu_path_p_vec1(i) + cumu_path_p_vec1(i-1);
end

%% linear array
path_idx = linspace(1, length(path_p_vec1), length(path_p_vec1));
path_idx = path_idx';

%% plot
fig = figure();
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
p1 = plot([path_idx(1), path_idx(end_idx1)], [const_c1, const_c1], ...
    'color',co(1, :), 'LineWidth', 2); hold on;
% cumulative
p2 = plot(path_idx(1:delta1:end_idx1), cumu_path_p_vec1(1:delta1:end_idx1) , ...
    'color',co(2, :), 'linestyle', '-', 'LineWidth', 2); hold on;
% placeholder for legend
p3 = plot(NaN, NaN, 'LineStyle', ':', 'LineWidth', 2, 'color', co(3, :));
%% conc
set(gca,'GridLineStyle','--');
xlabel('#Path', 'FontSize', 20);
ylabel('Concentration', 'FontSize', 20);
ylim([10^(1.0*log10(cumu_path_p_vec1(1))), 10^(0.998*(log10(const_c1)))]);
% ylim([0.0, 10^(0.9995*(log10(const_c1)))]);

yyaxis right
error_data1 = (const_c1 - cumu_path_p_vec1)./const_c1.*100;
semilogy(path_idx(1:delta1:end_idx1), error_data1(1:delta1:end_idx1), ... 
    'LineStyle', ':', 'LineWidth', 2, 'color', co(3, :)); hold on;
ylabel('%Relative Error', 'FontSize', 20);
% set(gca, 'ytick', []);

%% global settings
grid on;
xlim([1, end_idx1]);
leg_h = legend([p1; p2; p3],'EXACT','SOHR', 'Percentage Error');
set(leg_h, 'FontSize', 14, 'Box', 'off');
% [left, bottom, weight, height]
set(leg_h, 'Position', [0.275, 0.2, 0.5, 0.1])

a_x = gca;
t_x = a_x.XLim(1) + 0.465*(a_x.XLim(2) - a_x.XLim(1));
t_y = a_x.YLim(1) + 0.818*(a_x.YLim(2) - a_x.YLim(1));
text(t_x, t_y, spe_name);

%% Zoom in figure1
% create a new pair of axes inside current figure
z1_position = [.23 .45 .25 .25];
axes('position',z1_position);
box on; % put box around new pair of axes
plot([path_idx(1), path_idx(end_idx2)], [const_c1, const_c1], ...
    'color',co(1, :), 'LineWidth', 2); hold on;
% cumulative
semilogy(path_idx(1:delta2:end_idx2), cumu_path_p_vec1(1:delta2:end_idx2) , ...
    'color',co(2, :), 'linestyle', '-', 'LineWidth', 2); hold on;
axis tight;
grid on;
ylim([10^(1.0*log10(cumu_path_p_vec1(1))), 10^(0.999*(log10(const_c1)))]);
%% Error
yyaxis right
semilogy(path_idx(1:delta2:end_idx2), error_data1(1:delta2:end_idx2), ... 
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
z2_position = [.55 .45 .25 .25];
% create a new pair of axes inside current figure
axes('position',z2_position);

% plot a placeholder
plot(NaN, NaN);

set(gca, 'ytick', []);
yyaxis right
box on; % put box around new pair of axes
% cumulative
zf2_1 = semilogy(path_idx(1:delta3:end_idx3), path_p_vec1(1:delta3:end_idx3) , ...
    'color',co(2, :), 'linestyle', '-.', 'LineWidth', 2); hold on;
tick_pos = linspace(ceil(log10(path_p_vec1(end_idx3))),floor(log10(path_p_vec1(1))), ...
    1+floor(log10(path_p_vec1(1))) - ceil(log10(path_p_vec1(end_idx3))));
tick_pos = arrayfun(@(x) 10^x, tick_pos);
set(gca, 'ytick', tick_pos);
axis tight;
grid on;
% legend
leg_z2 = legend(zf2_1,['PATHWAY' newline 'PROBABILITY']);
set(leg_z2, 'FontSize', 8, 'Box', 'off');
% [left, bottom, weight, height]
set(leg_z2, 'Position', [z2_position(1) + z2_position(3)*0.3, ...
    z2_position(2) + z2_position(4)*0.5, 0.1, 0.1]);

%% save to file
figname = strcat('pathway_prob_concentration_',num2str(spe_idx), '_2.png');
print(fig, fullfile(file_dir, folder1, figname), '-r200', '-dpng');