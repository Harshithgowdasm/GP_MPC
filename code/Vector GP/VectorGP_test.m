
clear all;
close all;

%state 1 data
xtr1 = 10*sort(rand(50,1)); 
f = @(x) sin(x)+sqrt(x);
sd = 0.2;
ytr1 = f(xtr1) + randn(size(xtr1,1),1)*sd; 
xte1 = linspace(0,10,50)'; 

%state 2 data
xtr2 = 10*sort(rand(50,1)); 
f2 = @(x) 0.0002*exp(x);
ytr2 = f2(xtr2) + randn(size(xtr2,1),1)*sd; 
xte2 = linspace(0,10,50)';

%state 3 data
xtr3 = 10*sort(rand(50,1)); 
f3 = @(x) 0.2*x+0.5;
ytr3 = f3(xtr3) + randn(size(xtr3,1),1)*sd; 
xte3 = linspace(0,10,50)';

% to vector form
X = [xtr1;xtr2;xtr3];
Y = [ytr1;ytr2;ytr3];
Xte = [xte1;xte2;xte3];


 %predict
[ymu,ys2,nlZ] = VectorGP(X,Y,Xte,3,50);         % 3 states and each 50 data points



%plot
ym = f(xte1);  % actual ouputs
ym2 = f2(xte2);
ym3 = f3(xte3);
figure, hold on
plot(xte1,ym,'Color','k','LineWidth',2)
plot(xte2,ym2,'Color','b','LineWidth',2)
plot(xte3,ym3,'Color','g','LineWidth',2)
plot(xte1,ymu{1},'Color','r','LineWidth',2)
plot(xte2,ymu{2},'Color','c','LineWidth',2)
plot(xte3,ymu{3},'Color','y','LineWidth',2)
legend('y1 function-sin(x)+sqrt(x)','y2 function-0.0002*exp(x)','y3 function-0.2*x+0.5', 'y1 GP','y2 GP','y3 GP')


