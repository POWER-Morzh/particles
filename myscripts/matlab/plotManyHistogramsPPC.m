clear all;

figPNG = figure; % this figure is used to print plots in PNG files
imgPath = 'img/'; % folder in wich the program will save all PNG's
source_path = '../../../../'; % folder with vtk files
filenamelist = char('histogramL2Error', ...
	'histogramL2Error_Coordinate', ...
    'histogramMaxError', ...
    'histogramMaxError_Coordinate', ...
    'histogramMaxOffset', ...
    'histogramMaxOffset_Coordinate');
% filenamelist = char('histogramL2Error');

PPC_FOR_HISTOGRAM = [30 50 100];

for i = 1:size (filenamelist,1)
    for p_i = 1:length(PPC_FOR_HISTOGRAM)
        clf; % Very important, is it's not done then many plots will be on one graph

        filename = [deblank(filenamelist(i,:)) '-' num2str(PPC_FOR_HISTOGRAM(p_i))];
        histogramData = load([source_path filename]);
        boundary_data = histogramData(1,:); % get boundarys from file
        small_boundary = boundary_data(1);
        big_boundary = boundary_data(2);
        small_power = log10(small_boundary); % get the exponent
        big_power = log10(big_boundary);
        histogramData(1,:) = []; % get rid of first line(we don't need it anymore)
        % add one collumn and row to get good working pcolor...
        histogramData = [ histogramData zeros(size(histogramData,1),1) ];
        histogramData = [ histogramData; zeros(size(histogramData,2),1)' ];

        img = figure(figPNG);
        hold on;
        ylabel('iterations');
        xlabel('Logarithm to base 10');
        % Create axis
        x = [(small_power-1) small_power:(big_power-small_power)/(size(histogramData,2)-3):(big_power) (big_power+1)];
        y = 0:size(histogramData,1)-1;
        [X,Y] = meshgrid(x,y);
        pcolor(X,Y,histogramData);
        colorbar;
        % Save image
        imgName = [filename '.png'];
        print(img, '-dpng', [imgPath imgName]);
    end
end
