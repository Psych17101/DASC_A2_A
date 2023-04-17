close all

% disp('-------------------------------------')
% disp('Demo: Moments of a circle with r=10')
% disp('-------------------------------------')

% r=10; 
% th=(0:1:359)'; 
% xy=r*[cos(th*pi/180),sin(th*pi/180)];

b1 = 1.25/100;
b2 = 2.5/100;
b3 = 2/100;
t1 = 0.0005715;
t2 = 0.0005715;
t3 =  0.000762;


xy1 = [0 0];
xy2 = [b1/2 0];
xy3 = [b1/2, t1];
xy4 = [t2, t1];
xy5 = [t2, b2];
xy6 = [b3/2, b2];
xy7 = [b3/2, t3 + b2];
xy8 = [-b3/2,t3 + b2];
xy9 = [-b3/2,b2];
xy10 = [-t2,b2];
xy11 = [-t2,t1];
xy12 = [-b1/2,t1];
xy13 = [-b1/2,0];



xy = [xy1;xy2;xy3;xy4;xy5;xy6;xy7;xy8;xy9;xy10;xy11;xy12;xy13];

Polygon=PolygonMoments(xy,[],-3)
disp 'In this example, because of the 360 degree symmetry,'
disp 'any x-y axes are principal axes'

