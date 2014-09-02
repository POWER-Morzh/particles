clear all;

source_path = '../../../../';
filename = 'maxError-24.dat';

maxNorm = load([source_path filename]);

x = sort(maxNorm(:, 2)');
y = sort(maxNorm(:, 3)');

x = unique(x, 'R2012a');
y = unique(y, 'R2012a');

x = [x y];

x = unique(x, 'R2012a');
y = x;

[X Y] = meshgrid(1:10, 1:10);ZI = griddata([1 2 3 4],[1 2 3 6],[1 2 1 2],X,Y);mesh(ZI) ;

[X Y] = meshgrid(x, y);

Z = griddata(maxNorm(:,2)', maxNorm(:,3)', maxNorm(:,1)', X, Y);

mesh(X,Y,Z);


