clear all;

l2Norm = load('../coordinateL2ErrorNorm-9.dat');

figure;

x = sort(l2Norm(:, 3)');
y = sort(l2Norm(:, 4)');

x = unique(x, 'R2012a');
y = unique(y, 'R2012a');

x = [x y];

x = unique(x, 'R2012a');
y = x;

[X Y] = meshgrid(1:10, 1:10);ZI = griddata([1 2 3 4],[1 2 3 6],[1 2 1 2],X,Y);mesh(ZI) ;

[X Y] = meshgrid(x, y);

Z = griddata(l2Norm(:,3)', l2Norm(:,4)', l2Norm(:,1)', X, Y);

mesh(X,Y,Z);


