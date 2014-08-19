clear all;

histogramData = load('../../../../histogramL2ErrorOffset');
% get boundarys from file
boundary_data = histogramData(1,:);
small_boundary = boundary_data(1);
big_boundary = boundary_data(2);
% get the power
small_power = log10(small_boundary);
big_power = log10(big_boundary);
% get rid of first line(we don't need it anymore)
histogramData(1,:) = [];

% add one collumn and row to get good working pcolor...
histogramData = [ histogramData zeros(size(histogramData,1),1) ];
histogramData = [ histogramData; zeros(size(histogramData,2),1)' ];

figure;
hold on;

ylabel('iterations');
xlabel('Logarithm to base 10');



x = [(small_power-1) small_power:(big_power-small_power)/(size(histogramData,2)-3):(big_power) (big_power+1)];
y = 0:size(histogramData,1)-1;

[X,Y] = meshgrid(x,y);

pcolor(X,Y,histogramData);

%image(x,y,histogramData);