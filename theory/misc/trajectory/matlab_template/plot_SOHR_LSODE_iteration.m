function plot_SOHR_LSODE_iteration(speIndex)
% speIndex = 8;

addpath('./jsonlab-1.5');
%% Current file directory
file_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..');

%% dlsode time
filename = fullfile(file_dir, 'output', 'time_dlsode_fraction.csv');
delimiter = '';
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
dlsode_time = dataArray{:, 1};
clearvars filename delimiter formatSpec fileID dataArray ans;

%% dlsode concentration
filename = fullfile(file_dir, 'output', 'concentration_dlsode_fraction.csv');
delimiter = ',';
formatSpec = '%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
dlsode_concentration = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray ans;

%% SOHR time
filename = fullfile(file_dir, 'output', 'time_SOHR_fraction_all.csv');
delimiter = '';
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
SOHR_time = dataArray{:, 1};
clearvars filename delimiter formatSpec fileID dataArray ans;

%% iterationNumber
data = loadjson(fullfile(file_dir, 'input', 'setting.json'));
iterationNumber = data.SOHR_init.iterationNumber;

%% anchor time points
timeN1 = data.SOHR_init.timeN1;
tn = floor(length(SOHR_time)/timeN1);
anchor_time_point_index = (linspace(1, tn, tn).*timeN1);

%% SOHR concentration
SOHR_concentration = cell(iterationNumber, 1);
for i=1:iterationNumber
    filename = fullfile(file_dir, 'output', ['concentration_SOHR_fraction_all_', num2str(i), '.csv']);
    delimiter = ',';
    formatSpec = '%f%f%f%f%f%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
    fclose(fileID);
    SOHR_concentration{i, 1} = [dataArray{1:end-1}];
    clearvars filename delimiter formatSpec fileID dataArray ans; 
end

%% plot
% step1 = 10;
step1 = 1;
step2 = 1;
% [0-->O2, 1-->H2O, 2-->H2, 3-->H2O2, 4-->H, 5-->OH, 6-->HO2, 7-->O]
speLabel = {'$O_2$', '$H_2O$', '$H_2$', '$H_2O_2$', '$H$', '$OH$', '$HO_2$', '$O$'};
speName = {'O2', 'H2O', 'H2', 'H2O2', 'H', 'OH', 'HO2', 'O'};
% speIndex = 1;

colors = jet(iterationNumber + 1);
marker = {'hexagram', 'o', '*', '.', 'x', 'square', 'diamond', '^', 'v', '>', '<', 'pentagram', '+'};

fig = figure();
% Handler array
H = gobjects(iterationNumber + 2);
color_counter = 1;
marker_counter = 1;
initialGuessCounter = 1;
handler_counter = 1;

for itr = 1:step2:iterationNumber
    startIndex = 1;
    for endIndex = anchor_time_point_index
        % Initial Guess
        if itr == 1
            yDataPoint = SOHR_concentration{iterationNumber,1}(startIndex, speIndex) - dlsode_concentration(1,speIndex);
            iHandler = plot([SOHR_time(startIndex), SOHR_time(endIndex+1)], ...
                [yDataPoint, yDataPoint], ... 
                'color', colors(initialGuessCounter, :));
            hold on;
        end
        % iterations
        yData = SOHR_concentration{itr,1}(startIndex:step1:endIndex+1, speIndex) - dlsode_concentration(1,speIndex);
        if endIndex > anchor_time_point_index(1)
            yData(1) = SOHR_concentration{iterationNumber,1}(startIndex, speIndex) - dlsode_concentration(1,speIndex);
        end
        pHandler = plot(SOHR_time(startIndex:step1:endIndex+1), yData, ... 
            'color', colors(color_counter+1, :));
        hold on;
        % add legends, only for the first interval
        if endIndex == anchor_time_point_index(1)
            if itr == 1
                H(handler_counter) = iHandler;
                handler_counter = handler_counter + 1;
            end 
            H(handler_counter) = pHandler;
            handler_counter = handler_counter + 1;
        end          
        startIndex = endIndex + 1;
    end
    color_counter = color_counter + 1;
    
end

%% plot lsode
ratio = 1.0;
NExact = 10;

xData = dlsode_time(1:NExact:floor(length(dlsode_time)*ratio));
yData = dlsode_concentration(1:NExact:floor(length(dlsode_concentration(:, speIndex))*ratio),speIndex) - dlsode_concentration(1,speIndex);
H(handler_counter) = scatter(xData, yData, ...
    'marker', marker{1, marker_counter}, ...
    'MarkerFaceColor', colors(color_counter, :), ...
    'MarkerEdgeColor', colors(color_counter, :));
hold on;

%% configuration
xlim([dlsode_time(1), dlsode_time(floor(length(dlsode_time)*ratio))]);
grid on;

% xlabel('Time', 'FontSize', 20);
% ylabel('Concentration', 'FontSize', 20);
title(speLabel(speIndex), 'FontSize', 20, 'Interpreter','latex');

yt = get(gca, 'ytick');
yaxis = get(gca, 'yaxis');
% get the exponent of yticks
expVal = double(yaxis.Exponent);

if expVal ~= 0
    yt_label = strsplit(sprintf('%.1f,',yt/10^expVal), ',');
    set(gca, 'yticklabel', yt_label);  
end

pos = get(gca, 'Position');
x_offset = pos(3)*0.0;
y_offset = - pos(4)*0.015;
if dlsode_concentration(1,speIndex)>0
    append_str = ['+',num2str(dlsode_concentration(1,speIndex))];
else
    append_str='';
end

if expVal == 0
    annotation('textbox',[pos(1)+x_offset, pos(2)+pos(4)+y_offset, 0.2, 0.2],...
    'String', append_str,...
    'Interpreter','latex',...
    'VerticalAlignment','bottom',...
    'EdgeColor','none')
else
    annotation('textbox',[pos(1)+x_offset, pos(2)+pos(4)+y_offset, 0.2, 0.2],...
    'String',['$\times10^', '{', num2str(expVal), '}', append_str, '$'],...
    'Interpreter','latex',...
    'VerticalAlignment','bottom',...
    'EdgeColor','none')
end

% use lambda expression
legName = cell(1, 2+iterationNumber);
legName{1, end} = 'EXACT';
legName{1, 1} = 'Initial Guess';
for itrIndex = 1:iterationNumber
    legName{1, itrIndex+1} = ['n = ', num2str(itrIndex)];
end 

legHandler = legend(H(:,1) ,legName);
set(legHandler, 'FontSize', 14, 'Box', 'off');
if dlsode_concentration(1,speIndex) > 0
    set(legHandler, 'Location', 'SouthWest');
else
    set(legHandler, 'Location', 'NorthWest');
end
% h_pos=get(legHandler,'position');
% set(legHandler, 'position', [h_pos(1)*0.7, h_pos(2)*0.8, h_pos(3), h_pos(4)]);

%% save to file
fname = strjoin(strcat('Iterations_', speName(speIndex), '.png'));
print(fig, fullfile(file_dir, 'output', fname), '-r200', '-dpng');