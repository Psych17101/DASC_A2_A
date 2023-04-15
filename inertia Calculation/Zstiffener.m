close all

% disp('-------------------------------------')
% disp('Demo: Moments of a circle with r=10')
% disp('-------------------------------------')

% r=10; 
% th=(0:1:359)'; 
% xy=r*[cos(th*pi/180),sin(th*pi/180)];
b1 = 2/100;
b2 = 2/100;
b3 = 3/100;
t = 0.0762/100;

xy1 = [0 0];
xy2 = [b1-t, 0];
xy3 = [b1-t, t];
xy4 = [0, t];
xy5 = [0, b3 + t];
xy6 = [-b2, b3 + t];
xy7 = [-b2, b3];
xy8 = [-t,b3 ];
xy9 = [-t,0 ];

xy = [xy1;xy2;xy3;xy4;xy5;xy6;xy7;xy8;xy9];

Polygon=PolygonMoments(xy,[],-4)
disp 'In this example, because of the 360 degree symmetry,'
disp 'any x-y axes are principal axes'

