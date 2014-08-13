clear all;

figPNG = figure('units','normalized','outerposition',[0 0 1 1]); % this figure is used to print plots in PNG files

imgPath = '../img/'; % folder in wich the program will save all PNG
filePath = '../output/'; % folder with dat files

FACTOR = 2.^[0:0];
MANTISSA = 7:8;
FACTOR = [0.5 FACTOR];

filenamelist = char('coordinateGlobalL2ErrorNorm', ...
	'coordinate_global_max_error', ...
    'coordinateGlobalMaxL2ErrorNorm', ...
    'globalL2ErrorNorm', ...
    'global_max_error', ...
    'globalMaxL2ErrorNorm');
% filenamelist = char('coordinateGlobalL2ErrorNorm');

% for i = 1:size (filenamelist,1)
%     for f_i = 1:length(FACTOR)
%         filename = [deblank(filenamelist(i,:))];
%         hfigSurf(f_i,i) = figure( 'units','normalized','outerposition',[0 0 1 1],...
%             'Name',[filename '-FACTOR-' num2str(FACTOR(f_i))] );% Creation of figures for different times
%     end
% end


for i = 1:size (filenamelist,1)
    for f_i = 1:length(FACTOR)
        savefilename = [ deblank(filenamelist(i,:)) '-F' num2str(FACTOR(f_i)) ];
        clf; % Very important, is it's not done then many plots will be on one graph
        for m_i = 1:length(MANTISSA)
            filename = [ savefilename '-M' num2str(MANTISSA(m_i)) ];

            globalNorm = load([filePath filename '.dat']);
            
            
            %figure(hfigSurf(f_i,i)); % Create figure which will be divided in 4 parts
            img = figure(figPNG);
                subplot(2,2,m_i);
                hold on;
                xlabel(['MANTISSA='  num2str(MANTISSA(m_i))]);
                plot( globalNorm(:,1) );
                hold off;
                
            imgName = [savefilename '.png'];
            print(img, '-dpng', [imgPath imgName]);

            
        end
    end
end

