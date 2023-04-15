close all

% disp('-------------------------------------')
% disp('Demo: Moments of a circle with r=10')
% disp('-------------------------------------')

% r=10; 
% th=(0:1:359)'; 
% xy=r*[cos(th*pi/180),sin(th*pi/180)];

b1 = 0.15;
t = 0.9/1000;
b2 = 0.4-2*t;
b3 = 0.15;

xy1 = [0 0];
xy2 = [b1, 0];
xy3 = [b1, t];
xy4 = [t, t];
xy5 = [t, b2];
xy6 = [b1, b2];
xy7 = [b1, b2 + t];
xy8 = [0,b2+t];

xy = [xy1;xy2;xy3;xy4;xy5;xy6;xy7;xy8];

Polygon=PolygonMoments(xy,[],0.001)
disp 'In this example, because of the 360 degree symmetry,'
disp 'any x-y axes are principal axes'

