clear all;
figPNG = figure('units','normalized','outerposition',[0 0 1 1]); % this figure is used to print plots in PNG files
imgPath = 'pic/'; % folder in wich the program will save all PNG's
source_path = '../../../../output/'; % folder with vtk files
copy_in_source_path = 'img/'; % Where copy all files

filenamelist = char('histogramL2Error', ...
	'histogramL2Error_Coordinate', ...
    'histogramMaxError', ...
    'histogramMaxError_Coordinate', ...
    'histogramMaxOffset', ...
    'histogramMaxOffset_Coordinate');
% filenamelist = char('histogramL2Error');

PPC_FOR_HISTOGRAM = [30 50 100];
% Read first file from list to know which foldername to give
filename = [deblank(filenamelist(1,:)) '-' num2str(PPC_FOR_HISTOGRAM(1))];
histogramData = load([source_path filename]);
file_data = histogramData(1,:); % get first line
maximal_velocity = file_data(4);
foldername = [copy_in_source_path '250000-0.1-' sprintf('%.2f',maximal_velocity) '/'];
system(['mkdir -p ' foldername]); % Invoke system command(mkdir)
system(['mkdir -p ' foldername imgPath]);

for i = 1:size (filenamelist,1)
    savefilename = deblank(filenamelist(i,:));
    clf; % Very important, is it's not done then many plots will be on one graph
    for p_i = 1:length(PPC_FOR_HISTOGRAM)
        filename = [deblank(filenamelist(i,:)) '-' num2str(PPC_FOR_HISTOGRAM(p_i))];
        % Copy "vtk" data to save it for next uses
        system(['cp ' source_path filename ' ' foldername filename]);
        
        histogramData = load([source_path filename]);
        file_data = histogramData(1,:); % get boundarys from file
        small_boundary = file_data(1);
        big_boundary = file_data(2);
        small_power = log10(small_boundary); % get the exponent
        big_power = log10(big_boundary);
        histogramData(1,:) = []; % get rid of first line(we don't need it anymore)
        % add one collumn and row to get good working pcolor...
        histogramData = [ histogramData zeros(size(histogramData,1),1) ];
        histogramData = [ histogramData; zeros(size(histogramData,2),1)' ];

        img = figure(figPNG);
            subplot(2,2,p_i);
            hold on;
            ylabel('iterations');
            xlabel(['PPC='  num2str(PPC_FOR_HISTOGRAM(p_i))]);
            % Create axis
            x = [(small_power-1) small_power:(big_power-small_power)/(size(histogramData,2)-3):(big_power) (big_power+1)];
            y = 0:size(histogramData,1)-1;
            [X,Y] = meshgrid(x,y);
            pcolor(X,Y,histogramData);
            colorbar;
        % Save image
        imgName = [savefilename '.png'];
        print(img, '-dpng', [foldername imgPath imgName]);
    end
end
