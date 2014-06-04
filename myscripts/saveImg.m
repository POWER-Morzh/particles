clear all;

figPNG = figure; % this figure is used to print plots in PNG files

imgPath = '../img/'; % folder in wich the program will save all PNG
filePath = '../output/'; % folder with dat files

FACTOR = 2.^[0:1];
MANTISSA = 7:8;

filenamelist = char('globalL2ErrorNorm', 'global_max_error', 'globalMaxL2ErrorNorm', 'coordinateGlobalL2ErrorNorm', 'coordinate_global_max_error', 'coordinateGlobalMaxL2ErrorNorm');
%filenamelist = char('globalL2ErrorNorm','global_max_error');

for f = FACTOR
    for m = MANTISSA
        for i = 1:size (filenamelist,2)
            img = figure(figPNG);

            filename = [deblank(filenamelist(i,:)) '-F' num2str(f) '-M' num2str(m)];

            globalNorm = load([filePath filename '.dat']);
            plot( globalNorm(:,1) );

            imgName = [filename '.png'];
            print(img, '-dpng', [imgPath imgName]);
        end
    end
end