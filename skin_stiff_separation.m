%% skin-stiffener separation S 9.2.2
clear; format;
close;
format short g;

% Define material properties
E_x = 62.046*10^9;
E_y = 62.046*10^9;
E1   = 165*10^9;% Pa - direction modulus of the lamina
E2   = 8*10^9  ; % Pa
nu12 = 0.05 ;
nu21 = nu12*E1/E2;
G12  = 6*10^9  ; % Pa
X_T = 1517*10^6;
X_C = 1379*10^6;
Y_T = 1450*10^6;
Y_C = 1379*10^6;
S = 99*10^6;


% Panel dimensions
a = 1; % [m] y direction 
bp = 0.48; % [m] x direction (side where force is applied)
n_S=5; %number of stiffeners
ds = bp/n_S; 
b_stiff=5.5e-2;
%% the interlaminar stresses die out after 10-15 flange thicknesses away from
%the flanged edge p.249, 

t_ply=0.1905e-3/2; %for the convention of (0/90)=1 thickness
t=t_ply*8; %thickness of the skin
death=15*t;
next=(bp-b_stiff*n_S)/(n_S-1);
dontworry=death<next; %if this is true, it dies out before the next stiffener

%if this is alright, the interaction between the stiffeners can be
%neglected, nr.1 of the guidelines from page 249 is achieved
%%  

















