close all

% disp('-------------------------------------')
% disp('Demo: Moments of a circle with r=10')
% disp('-------------------------------------')

% r=10; 
% th=(0:1:359)'; 
% xy=r*[cos(th*pi/180),sin(th*pi/180)];

b= 2/100;

t= 0.0762/100;
l1 = 0.25/100;
l2 = 1.5/100;
theta = acos(l1/l2);
h= l2*sin(theta);
xy1 = [0 0];
xy2 = [b/2 0];
xy3 = [b/2 + l1, -h];
xy4 = [b/2 + l1 + l2, -h];
xy5 = [b/2 + l1 + l2, -h-t];
xy6 = [b/2 + l1, -h-t];
xy7 = [b/2 -t];
xy8 = [-b/2 -t];
xy9 = [-b/2 - l1, -h-t];
xy10 = [-b/2 - l1 - l2, -h-t];
xy11 = [-b/2 - l1 - l2, -h];
xy12 = [-b/2 - l1, -h];
xy13 = [-b/2 0];

xy = [xy1;xy2;xy3;xy4;xy5;xy6;xy7;xy8;xy9;xy10;xy11;xy12;xy13];

Polygon=PolygonMoments(xy,[],-4)
disp 'In this example, because of the 360 degree symmetry,'
disp 'any x-y axes are principal axes'



