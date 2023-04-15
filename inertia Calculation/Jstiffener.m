close all

% disp('-------------------------------------')
% disp('Demo: Moments of a circle with r=10')
% disp('-------------------------------------')

% r=10; 
% th=(0:1:359)'; 
% xy=r*[cos(th*pi/180),sin(th*pi/180)];

b3 = 3/100;
b2 = 3/100;
b1 = 1/100;
t = 0.0762/100;

xy1 = [0 0];
xy2 = [b3/2 0];
xy3 = [b3/2, t];
xy4 = [t/2, t];
xy5 = [t/2, b2 + t];
xy6 = [-b1 - t/2, b2 + t];
xy7 = [-b1 - t/2, b2];
xy8 = [-t/2 , b2];
xy9 = [-t/2, t];
xy10 = [-b3/2, t];
xy11 = [-b3/2 , 0];


xy = [xy1;xy2;xy3;xy4;xy5;xy6;xy7;xy8;xy9;xy10;xy11];

Polygon=PolygonMoments(xy,[],-4)
disp 'In this example, because of the 360 degree symmetry,'
disp 'any x-y axes are principal axes'

