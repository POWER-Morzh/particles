clear all;
figure;
hold on;
xlabel('timestep');
globalNorm = load('../globalL2ErrorNorm');


plot( globalNorm(:,1) );