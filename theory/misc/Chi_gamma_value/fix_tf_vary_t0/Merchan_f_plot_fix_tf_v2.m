%% global settings
file_dir = fullfile(fileparts(mfilename('fullpath')));

tau = 0.777660157519;
end_t = '0.9';
end_t2 = 0.9;

%% time
fn_time = fullfile(file_dir, ['pathway_time_candidate_S10_HA4_', end_t ,'.csv']);
delimiter = ',';
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(fn_time,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
time_mat = [dataArray{1:end-1}];
time_vec = time_mat(1, :);
clearvars fn_time delimiter formatSpec fileID dataArray ans time_mat;

%% pathwap probability
fn_pp = fullfile(file_dir, ['pathway_prob_S10_HA4_', end_t ,'.csv']);
delimiter = ',';
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(fn_pp,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
pp_mat = [dataArray{1:end-1}];
clearvars fn_pp delimiter formatSpec fileID dataArray ans;

%% total top 3 path
top_3_idx = [1, 2, 3];
for idx=1:length(top_3_idx)
    if idx == 1
        top_3_vec = pp_mat(idx, :);
    else
        top_3_vec = (top_3_vec + pp_mat(idx, :));
    end
end

%% total top 1000 path
top_4_to_1000_idx = linspace(4, 1000, 1000 - 3);
for idx=1:length(top_4_to_1000_idx)
    true_idx = top_4_to_1000_idx(idx);
    if idx == 1
        top_4_to_1000_vec = pp_mat(true_idx, :);
    else
        top_4_to_1000_vec = (top_4_to_1000_vec + pp_mat(true_idx, :));
    end
end

%% top3_over_top1000
top3_over_top4to1000_vec = top_3_vec ./ top_4_to_1000_vec;


%% plot
%% plot
fig = figure();
% https://www.mathworks.com/help/matlab/graphics_transition/why-are-plot-lines-different-colors.html
% https://www.mathworks.com/help/matlab/creating_plots/customize-graph-with-two-y-axes.html
co = [    0    0.4470    0.7410 % 1th plot
    0.8500    0.3250    0.0980 % 2nd plot
    0.9290    0.6940    0.1250 % 3rd plot
    0.4940    0.1840    0.5560 % 4th plot
%     0.4660    0.6740    0.1880 % 5th plot
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
p2= plot(time_vec * tau, pp_mat(2, :), 'LineWidth', 2, 'marker', '>'); hold on;
p3 = plot(time_vec * tau, pp_mat(3, :), 'LineWidth', 2, 'marker', 's'); hold on;

p4 = plot(time_vec * tau, top_3_vec, 'LineWidth', 2, 'marker', 'd'); hold on;

p5 = plot(time_vec * tau, top_4_to_1000_vec, 'LineWidth', 2, 'marker', 'o'); hold on;
% p5 = plot(nan, nan, 'LineWidth', 2, 'marker', 'd');

%% settings
set(gca,'GridLineStyle','--');
xlabel('Time (seconds)', 'FontSize', 20);
ylabel('f', 'FontSize', 20);
% ylim([0.15, 1.8]);

% %% ratio
% yyaxis right
% delta_n = 1000;
% plot(time_vec * tau, top3_over_top4to1000_vec, 'LineWidth', 2, 'color', 'r', 'marker', 'd'); hold on;
% ylabel('Ratio', 'FontSize', 20);
% % set(gca, 'ytick', []);

%% global settings
grid on;
xlim([0, tau*end_t2]);
leg_h = legend([p1; p2; p3; p4; p5], 'p1', 'p2', 'p3', 'SUM(p1,p2,p3)', 'SUM(p4,...,p1000)');
set(leg_h, 'FontSize', 12, 'Box', 'off');
% set(leg_h, 'Location', 'South')

%% save to file
figname = strcat('Merchant_f_fix_tf_vary_t0_', end_t, '_v2.png');
print(fig, fullfile(file_dir, figname), '-r200', '-dpng');



