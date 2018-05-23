%% global settings

%% Current file directory
file_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', '..', '..', '..', '..', '..', '..', 'SOHR_DATA');
pic_dir = fullfile(fileparts(mfilename('fullpath')));

%% import time
fn_time = fullfile(file_dir, 'output', 'time_dlsode_M.csv');
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
fn_temp = fullfile(file_dir, 'output', 'temperature_dlsode_M.csv');
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
fn_R = fullfile(file_dir, 'output', 'reaction_rate_dlsode_M.csv');
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

%% plot
fig = figure();
% https://www.mathworks.com/help/matlab/graphics_transition/why-are-plot-lines-different-colors.html
% https://www.mathworks.com/help/matlab/creating_plots/customize-graph-with-two-y-axes.html
co = [    0    0.4470    0.7410 % 1th plot
    0.8500    0.3250    0.0980 % 2nd plot
%     0.9290    0.6940    0.1250 % 3rd plot
%     0.4940    0.1840    0.5560 % 4th plot
%     0.4660    0.6740    0.1880 % 5th plot
%     0.3010    0.7450    0.9330 % 6th plot
%     0.6350    0.0780    0.1840 % 7th plot
%     0   0   1 % placeholder
%     0   0.5   0 % placeholder
%     0   0.75   0.75 % placeholder
%     0.7500         0    0.7500 % placeholder
%     0.7500    0.7500         0 % placeholder
%     0.2500    0.2500    0.2500 % placeholder
    1   0   0 % bl
    ]; 
set(fig,'defaultAxesColorOrder',co)
%%
tau = 0.777660157519;
end_t = 0.9;

% data preparation alpha
R_idx_sink1 = [1, 3, 4, 7, 21, 22, 29, 31, 32, 37, 39, 43, 44, 51, 52, 61, 62, 67, 79, 89, 90, 104, 111, 121, 122, 149, 152, 155, 161, 163, 164, 179, 183, 185, 186, 198, 209, 210, 212, 228, 236, 241, 243, 262, 267, 282, 290, 311, 312, 339, 361, 367, 369, 373, 377, 379, 383, 388, 411, 412, 424, 431, 447, 453, 461, 462, 464, 466, 471, 482, 501, 514, 517, 522, 524, 526, 533, 542, 544, 546, 567, 569, 571, 582, 589, 594, 599, 623, 628, 653, 654, 680, 685, 715, 717, 733, 735, 736, 738, 786, 815, 817, 819, 853, 862, 867, 884, 889, 898, 903, 907, 909, 1002, 1007, 1014, 1016, 1018, 1020, 1022, 1024, 1026, 1028, 1035, 1039, 1077, 1085, 1091, 1095, 1103, 1109, 1113, 1115, 1119, 1123, 1131, 1133, 1135, 1151, 1153, 1155, 1163, 1167, 1173, 1175, 1177, 1187, 1189, 1191, 1199, 1201, 1203, 1211, 1213, 1215, 1217, 1219, 1221];
for idx=1:length(R_idx_sink1)
    r_idx = R_idx_sink1(idx) + 1;
    if idx == 1
        R_total1 = reaction_R_mat(:, r_idx);
    else
        R_total1 = (R_total1 + reaction_R_mat(:, r_idx));
    end
end
% R_idx_source1 = [0, 2, 5, 6, 20, 23, 28, 30, 33, 36, 38, 42, 45, 50, 53, 60, 63, 66, 78, 88, 91, 105, 110, 120, 123, 148, 153, 154, 160, 162, 165, 178, 182, 184, 187, 199, 208, 211, 213, 229, 237, 240, 242, 263, 266, 283, 291, 310, 313, 338, 360, 366, 368, 372, 376, 378, 382, 389, 410, 413, 425, 430, 446, 452, 460, 463, 465, 467, 470, 483, 500, 515, 516, 523, 525, 527, 532, 543, 545, 547, 566, 568, 570, 583, 588, 595, 598, 622, 629, 652, 655, 681, 684, 714, 716, 732, 734, 737, 739, 787, 814, 816, 818, 852, 863, 866, 885, 888, 899, 902, 906, 908, 1003, 1006, 1015, 1017, 1019, 1021, 1023, 1025, 1027, 1029, 1034, 1038, 1076, 1084, 1090, 1094, 1102, 1108, 1112, 1114, 1118, 1122, 1130, 1132, 1134, 1150, 1152, 1154, 1162, 1166, 1172, 1174, 1176, 1186, 1188, 1190, 1198, 1200, 1202, 1210, 1212, 1214, 1216, 1218, 1220];
% for idx=1:length(R_idx_source1)
%     r_idx = R_idx_source1(idx) + 1;
%     R_total1 = (R_total1 - reaction_R_mat(:, r_idx));
% end

r_sink1 = [736];
for idx=1:length(r_sink1)
    r_idx = r_sink1(idx)+1;
    if idx==1
        data_y1 = reaction_R_mat(:, r_idx);
    else
        data_y1 = (data_y1 + reaction_R_mat(:, r_idx));
    end
end
r_source1 = [737];
for idx=1:length(r_source1)
    r_idx = r_source1(idx)+1;
    data_y1 = (data_y1 - reaction_R_mat(:,  r_idx));
end

h1 = plot(time_vec, data_y1 ./ (R_total1), ...
    'LineWidth', 2); hold on;

% data preparation alpha
R_idx_sink2 = [706, 710, 856, 914, 916, 922, 926, 930, 934, 938, 942, 946, 952, 956, 960, 964, 968, 970, 982, 984, 986, 988, 990, 1069, 1078, 1080, 1082, 1084];
%% chattering reactions
% chattering_R_idx = [1068, 1069];
chattering_R_idx = [];
R_idx_sink2 = setdiff(R_idx_sink2, chattering_R_idx);
for idx=1:length(R_idx_sink2)
    r_idx = R_idx_sink2(idx) + 1;
    if idx == 1
        R_total2 = reaction_R_mat(:, r_idx);
    else
        R_total2 = (R_total2 + reaction_R_mat(:, r_idx));
    end
end
% R_idx_source2 = [707, 711, 857, 915, 917, 923, 927, 931, 935, 939, 943, 947, 953, 957, 961, 965, 969, 971, 983, 985, 987, 989, 991, 1068, 1079, 1081, 1083, 1085];
% R_idx_source2 = setdiff(R_idx_source2, chattering_R_idx);
% for idx=1:length(R_idx_source2)
%     r_idx = R_idx_source2(idx) + 1;
%     R_total2 = (R_total2 - reaction_R_mat(:, r_idx));
% end

% r_sink2 = [1080, 1078];
r_sink2 = [1080];
for idx=1:length(r_sink2)
    r_idx = r_sink2(idx)+1;
    if idx==1
        data_y2 = reaction_R_mat(:, r_idx);
    else
        data_y2 = (data_y2 + reaction_R_mat(:, r_idx));
    end
end
% r_source2 = [1081, 1079];
r_source2 = [1081];
for idx=1:length(r_source2)
    r_idx = r_source2(idx)+1;
    data_y2 = (data_y2 - reaction_R_mat(:,  r_idx));
end

data_plot2 = data_y2 ./ (R_total2);
h2 = plot(time_vec, data_plot2, ...
    'LineWidth', 2); hold on;

%% settings
set(gca,'GridLineStyle','--');
xlabel('$t$ (seconds)', 'FontSize', 20, 'Interpreter', 'Latex');
ylabel('Fraction', 'FontSize', 20);
% ylim([10^-3.5, 10^0]);

%% temp
yyaxis right
delta_n = 1000;
plot(time_vec, temp_vec, 'LineWidth', 2, 'color', 'r'); hold on;
pt = scatter(time_vec(1:delta_n:end), temp_vec(1:delta_n:end), 'MarkerEdgeColor', 'r');
ylabel('T (K)', 'FontSize', 20);
% set(gca, 'ytick', []);

%% global settings
grid on;
xlim([0, tau*end_t]);

leg_h = legend('\alpha', '\beta');
set(leg_h, 'FontSize', 12, 'Box', 'off');
set(leg_h, 'Location', 'South')


%% save to file
figname = strcat('merchant_alpha_beta', '.png');
print(fig, fullfile(pic_dir, figname), '-r200', '-dpng');