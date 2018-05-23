% function subplot_SOHR_LSODE_iteration(speIndex1, speIndex2, speIndex3, speIndex4)
speIndex1=2;
speIndex2=4;
speIndex3=6;
speIndex4=8;

speIndex = [speIndex1, speIndex2, speIndex3, speIndex4];

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
maxIterationNumber = data.SOHR_init.iterationNumber;
iterationNumber = 1:2:maxIterationNumber;

%% anchor time points
timeN1 = data.SOHR_init.timeN1;
tn = floor(length(SOHR_time)/timeN1);
anchor_time_point_index = (linspace(1, tn, tn).*timeN1);

%% SOHR concentration
SOHR_concentration = cell(length(iterationNumber), 1);
for i=1:length(iterationNumber)
    filename = fullfile(file_dir, 'output', ['concentration_SOHR_fraction_all_', num2str(iterationNumber(i)), '.csv']);
    delimiter = ',';
    formatSpec = '%f%f%f%f%f%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
    fclose(fileID);
    SOHR_concentration{i, 1} = [dataArray{1:end-1}];
    clearvars filename delimiter formatSpec fileID dataArray ans; 
end

%% settings
step1 = 1;
step2 = 1;
% [0-->O2, 1-->H2O, 2-->H2, 3-->H2O2, 4-->H, 5-->OH, 6-->HO2, 7-->O]
speLabel = {'$O_2$', '$H_2O$', '$H_2$', '$H_2O_2$', '$H$', '$OH$', '$HO_2$', '$O$'};
speName = {'O2', 'H2O', 'H2', 'H2O2', 'H', 'OH', 'HO2', 'O'};
% speIndex = 1;

colors = prism(length(iterationNumber) + 2);
marker = {'o', '*', 'hexagram', '.', 'x', 'square', 'diamond', '^', 'v', '>', '<', 'pentagram', '+'};

fig = figure();

for i = 1:length(speIndex)
    subplot(2,2,i);
    % Handler array
    H = gobjects(length(iterationNumber) + 2);
    color_counter = 1;
    marker_counter = 1;
    initialGuessCounter = 1;
    handler_counter = 1;

    for itr = 1:step2:length(iterationNumber)
        startIndex = 1;
        for endIndex = anchor_time_point_index
            % Initial Guess
            if itr == 1
                yDataPoint = SOHR_concentration{length(iterationNumber),1}(startIndex, speIndex(i));
                iHandler = plot([SOHR_time(startIndex), SOHR_time(endIndex+1)], ...
                    [yDataPoint, yDataPoint], ... 
                    'LineWidth', 1.5, ...
                    'color', 'k');
                hold on;
            end
            % iterations
            yData = SOHR_concentration{itr,1}(startIndex:step1:endIndex+1, speIndex(i));
            if endIndex > anchor_time_point_index(1)
                yData(1) = SOHR_concentration{length(iterationNumber),1}(startIndex, speIndex(i));
            end
            if colors(length(colors)+1-color_counter, :) == [1 1 0]
                color_counter = color_counter + 1;
            end
            pHandler = plot(SOHR_time(startIndex:step1:endIndex+1), yData, ... 
                'LineWidth', 1.5, ...
                'color', colors(length(colors)+1-color_counter, :));
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
    yData = dlsode_concentration(1:NExact:floor(length(dlsode_concentration(:, speIndex(i)))*ratio),speIndex(i));
    H(handler_counter) = scatter(xData, yData, ...
        15.0, ...
        'marker', marker{1, marker_counter}, ...
        'LineWidth',1.0, ...
        'MarkerEdgeColor', colors(length(colors)+1-color_counter, :));
    hold on;

    %% configuration
    xlim([dlsode_time(1), dlsode_time(floor(length(dlsode_time)*ratio))]);
    grid on;

    if i == 3 || i == 4
        xlabel('Time', 'FontSize', 14);
    end
    
    pos = get(gca, 'Position');
    x_offset = pos(3)*0.45;
    y_offset = - pos(4)*0.2;  
    
    annotation('textbox',[pos(1)+x_offset, pos(2)+pos(4)+y_offset, 0.2, 0.2],...
    'String', speLabel(speIndex(i)),...
    'Interpreter','latex',...
    'VerticalAlignment','bottom',...
    'EdgeColor','none')

    if i == 1
        % use lambda expression
        legName = cell(1, 2+length(iterationNumber));
        legName{1, end} = 'EXACT';
        legName{1, 1} = 'Initial Guess';
        for itrIndex = 1:length(iterationNumber)
            legName{1, itrIndex+1} = ['n = ', num2str(iterationNumber(itrIndex))];
        end 

        legHandler = legend(H(:,1) ,legName);
        set(legHandler, 'FontSize', 10, 'Box', 'off');
        set(legHandler, 'Location', 'West');

    end
end


%% save to file
fname = strjoin(strcat('Iterations_', speName(speIndex1), speName(speIndex2), speName(speIndex3), speName(speIndex4), '.png'));
print(fig, fullfile(file_dir, 'output', fname), '-r200', '-dpng');
