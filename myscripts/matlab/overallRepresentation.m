clear all;

source_path = 'img/';
folder_name_list = char( ...
    '250000-0.1-0.07-normal', ...
    '250000-0.1-0.07-shifted', ...
    '250000-0.1-0.07-simple', ...
    '250000-0.1-0.71-normal', ...
    '250000-0.1-0.71-shifted', ...
    '250000-0.1-0.71-simple');
PPC_FOR_HISTOGRAM = [30 50 100];
filenamelist = char('histogramMaxError');

result_matrix = [];

for i=1:size(folder_name_list,1)
    for ppc_i = 1:length(PPC_FOR_HISTOGRAM)
        filename = [deblank(filenamelist(1,:)) '-' num2str(PPC_FOR_HISTOGRAM(ppc_i))];
        histogram_data = load([source_path deblank(folder_name_list(i,:)) '/' filename]);
        file_data = histogram_data(1,:); % get first line
        small_power = log10(file_data(1)); % get the exponent
        big_power = log10(file_data(2));
        histogram_data(1,:) = []; % get rid of first line(we don't need it anymore)

        % number of cells summed up by iteretion and range
        total_number_of_cells = sum(histogram_data(:));
        sum_of_rows = sum(histogram_data, 1);
        result_column = [sum_of_rows'; total_number_of_cells];
        result_matrix = [result_matrix result_column];
    end
end
% Output in file
csvfilename = 'Result.csv';
fid =  fopen([source_path csvfilename], 'w');
fprintf(fid, ' ,');
for i=1:size(folder_name_list,1)
    fprintf(fid,[deblank(folder_name_list(i,:)) ',,,']);
end
fprintf(fid, '\nppc,');
for i=1:size(folder_name_list,1)
    for j=1:length(PPC_FOR_HISTOGRAM)
        fprintf(fid,[num2str(PPC_FOR_HISTOGRAM(j)) ',']);
    end
end
fprintf(fid, '\n');
precisions = char(...
    '<1e-7', ...
    '<1e-6', ...
    '<1e-5', ...
    '<1e-4', ...
    '<1e-3', ...
    '>1e-3', ...
    'total number of cells');
for i=1:size(result_matrix,1)
    fprintf(fid,[deblank(precisions(i,:)) ',']);
    for j=1:size(result_matrix,2)
        fprintf(fid,[num2str(result_matrix(i,j)) ',']);
    end
    fprintf(fid, '\n');
end

fclose(fid);
