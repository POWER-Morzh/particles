clear all;

source_path = 'img/';
folder_name = '250000-0.1-0.07-simple/';
PPC_FOR_HISTOGRAM = [30 50 100];
% filenamelist = char('histogramL2Error', ...
% 	'histogramL2Error_Coordinate', ...
%     'histogramMaxError', ...
%     'histogramMaxError_Coordinate', ...
%     'histogramMaxOffset', ...
%     'histogramMaxOffset_Coordinate');
filenamelist = char('histogramMaxError');


filename = [deblank(filenamelist(1,:)) '-' num2str(PPC_FOR_HISTOGRAM(1))];
histogram_data = load([source_path folder_name filename]);
file_data = histogram_data(1,:); % get first line
small_power = log10(file_data(1)); % get the exponent
big_power = log10(file_data(2));
histogram_data(1,:) = []; % get rid of first line(we don't need it anymore)

% number of cells summed up by iteretion and range
total_number_of_cells = sum(histogram_data(:));
sum_of_rows = sum(histogram_data, 1);
normalized_sum_of_rows = sum_of_rows ./ total_number_of_cells;
normalized_sum_of_rows_for_plot = [normalized_sum_of_rows; zeros(1,length(normalized_sum_of_rows))]; 
% Create mesh
x = 0:1;
y = [(small_power-1) small_power:(big_power-small_power)/(size(histogram_data,2)-3):(big_power) (big_power+1)];
[X,Y] = meshgrid(x,y);
pcolor(X,Y,normalized_sum_of_rows_for_plot');
