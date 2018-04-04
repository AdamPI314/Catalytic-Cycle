%% global settings
file_dir = fullfile(fileparts(mfilename('fullpath')));

spe_idx = '60';
tau = 0.777660157519;
end_t = '0.9';
end_t2 = 0.9;

%% time
fn_time = fullfile(file_dir, ['pathway_time_candidate_S',spe_idx ,'_HA4_', end_t ,'.csv']);
delimiter = ',';
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(fn_time,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
time_mat = [dataArray{1:end-1}];
time_vec = time_mat(1, :);
clearvars fn_time delimiter formatSpec fileID dataArray ans time_mat;

%% pathwap probability
fn_pp = fullfile(file_dir, ['pathway_prob_S',spe_idx ,'_HA4_', end_t ,'.csv']);
delimiter = ',';
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(fn_pp,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
pp_mat = [dataArray{1:end-1}];
clearvars fn_pp delimiter formatSpec fileID dataArray ans;

%% plot
%% plot
fig = figure();
% https://www.mathworks.com/help/matlab/graphics_transition/why-are-plot-lines-different-colors.html
% https://www.mathworks.com/help/matlab/creating_plots/customize-graph-with-two-y-axes.html
co = [    0    0.4470    0.7410 % 1th plot
    0.8500    0.3250    0.0980 % 2nd plot
    0.9290    0.6940    0.1250 % 3rd plot
    0.4940    0.1840    0.5560 % 4th plot
    0.4660    0.6740    0.1880 % 5th plot
%     0.3010    0.7450    0.9330 % 6th plot
%     0.6350    0.0780    0.1840 % 7th plot
%     0   0   1 % placeholder
%     0   0.5   0 % placeholder
%     0   0.75   0.75 % placeholder
%     0.7500         0    0.7500 % placeholder
%     0.7500    0.7500         0 % placeholder
%     0.2500    0.2500    0.2500 % placeholder
    1   0   0 % placeholder
    1   0   0 % bl
    ]; 
set(fig,'defaultAxesColorOrder',co)

%% plot
p1 = plot(time_vec * tau, pp_mat(1, :), 'LineWidth', 2, 'marker', '<'); hold on;

%% settings
set(gca,'GridLineStyle','--');
xlabel('$t_f$ (seconds)', 'Interpreter','latex', 'FontSize', 20);
ylabel('$f$', 'Interpreter','latex', 'FontSize', 20);
% ylim([0.15, 1.8]);

%% global settings
grid on;
xlim([0, tau*end_t2]);
leg_h = legend(p1, 'SUM(p1,p2,p3)');
set(leg_h, 'FontSize', 12, 'Box', 'off');
% set(leg_h, 'Location', 'South')
% leg_pos = [0.5 0.425 0.1 0.1];
% set(leg_h, 'Position', [leg_pos(1),leg_pos(2), leg_pos(3), leg_pos(4)]);

%% text
a_x = gca;
t_x = a_x.XLim(1) + 0.325*(a_x.XLim(2) - a_x.XLim(1));
t_y = a_x.YLim(1) + 0.858*(a_x.YLim(2) - a_x.YLim(1));
text(t_x, t_y, 'npropyl, $t_0$=0', 'Interpreter','latex', 'FontSize', 20);

%% save to file
figname = strcat('2d_Merchant_f_', end_t, '_S', spe_idx, '_v1.png');
print(fig, fullfile(file_dir, figname), '-r200', '-dpng');



