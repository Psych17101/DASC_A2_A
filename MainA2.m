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
a = 0.48; % [m] y direction 
b = 1; % [m] x direction (side where force is applied)
n_stiff=5; %number of stiffeners
ds = a/n_stiff; 
%% Initialisation of Composites Model
% ABD Skin/Panel
% Laminate definition (plies of equal thickness)
%thetadt_skin = [0 90 +45 -45 -45 +45 90 0];%Original 
thetadt_skin = [0 90 +45 -45 -45 +45 90 0];% New tested
[A_skin,B_skin,D_skin,ABD_skin,h_skin,Qbar] = ABD_matrixCal(thetadt_skin,E_x,E_y,nu12,G12);
% Calcuation of global strain for first ply failure

F = [-4e3/a;0;0];
strain_glo = A_skin\F;
% Calculation of global stresses
stress_glo = Qbar(:,:,1)*strain_glo; % global sigmaxx etc...
Nplies = length(thetadt_skin);
for l = 1:Nplies
            % Calculations of local Strains and Stresses
            [eps_loc] = strain_gtol(strain_glo,thetadt_skin(l));
            [sigma_loc] = stress_gtol(stress_glo,thetadt_skin(l));% ply i angle in radians, from bottom
            % Failure index with Maximum Stress criterion
            [FI_1(l),FI_2(l),FI_3(l)]= MaxStress(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,S);
end

% ABD hat Stiffeners
%thetadt_hat = [0 90 0 90 90 0 90 0];% original
thetadt_hat = [0 90 0 90 90 0 90 0];% new test
[A_hat,B_hat,D_hat,ABD_hat,h_hat,Qbar_hat] = ABD_matrixCal(thetadt_hat,E_x,E_y,nu12,G12);

F = [-3e3/a;0;0];
strain_glo = A_hat\F;
stress_glo = Qbar_hat(:,:,1)*strain_glo; % global sigmaxx etc...
Nplies = length(thetadt_skin);
for l = 1:Nplies
            % Calculations of local Strains and Stresses
            [eps_loc] = strain_gtol(strain_glo,thetadt_skin(l));
            [sigma_loc] = stress_gtol(stress_glo,thetadt_skin(l));% ply i angle in radians, from bottom
            % Failure index with Maximum Stress criterion
            [FI_1(l),FI_2(l),FI_3(l)]= MaxStress(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,S);
end

%Equivalent Stiffness and Bending
A_hatT = n_stiff*A_hat;
A_eq = A_hatT + A_skin;

inv_A_hatT = inv(A_hatT);
EA_hat1 = 1/(inv_A_hatT(1,1)*h_hat);

D_hatT = n_stiff*D_hat;
D_eq = D_hatT + D_skin;

inv_D_hatT = inv(D_hatT);
EI_hat1 = 1/inv(D_hat(1,1)*h_skin);

%% Skin calulations % Section 8.8
D_skin(1,1) = 659.7/10^3;
D_skin(1,2) = 466.9/10^3;
D_skin(2,2) = 659.7/10^3;
D_skin(3,3) = 494/10^3;
AR_bar = b/ds;
ds = 0.457/4;
AR_bar = 50.8/ds;
invA_skin = inv(A_skin);
invA_hat = inv(A_hat);
E_skin = 1/(invA_skin(1,1)*h_skin); % equ 8.6
E_hat = 1/(invA_hat(1,1)*h_hat); % equ 8.6 - membrane stiffness of only on  
N_cr = pi^2/b^2*(D_skin(1,1)+2*(D_skin(1,2)+2*D_skin(3,3))+D_skin(2,2))/(a*(1+A_skin(1,2)/A_skin(1,1)));
N_x = 7*10^3/0.48;
% equa 7.15
%a_eff = a/(2*(1+2*(1+A_skin(1,2)/A_skin(1,1))*(N_cr/N_x)*(A_skin(1,1)/(A_skin(1,1)+3*A_skin(2,2)))));
a_eff = 0.303*a;
m = 10; %m_optimal 
b = 0.508;
N_x_skin = (pi^2/b^2)*(D_skin(1,1)*m^2 + 2*(D_skin(1,2)+2*D_skin(3,3))*AR_bar^2 + D_skin(2,2)*AR_bar^4/m^2); %(9.20)
% We know that with 7kN that our skin fails 
%N/m
Force_buckling = N_x_skin*a; %N 

%% Properties of Stiffeners
%hat
Area_hat = 4.191e-5;%Cross-sectional area of hat stiffener [m^2]
Interia_hat =  2.0793e-09;
members_hat = 5;
%I
Area_I = 8.2078e-05;%Cross-sectional area of hat stiffener [m^2]
Interia_I =  2.0394e-08;
members_I = 3;
%J
Area_J = 5.334e-05;%Cross-sectional area of hat stiffener [m^2]
Interia_K =  1.4432e-08;
members_J = 3;

EI_hat = members_hat*E_hat*Interia_hat;
Area_skin = a_eff*h_skin;
EA_hat = members_hat*E_hat*Area_hat; 
EA_skin = E_skin*Area_skin;
EA_total = EA_hat + EA_skin; % EA of one section stiffeners and skin

% Ratio of force on each member
rF_m1 = EA_hat/EA_total;
rF_m2 = EA_skin/EA_total;

%F_total_hat = rF_m1*Loadbeyond_buckling;
%F_total_skin = rF_m2*Load_buckling + Loadbeyond_buckling;
% 
% Applied_sigma_hat = F_total_skin/Area_hat;
% Applied_sigma_skin = F_total_skin/Area_skin;

lambda = (A_skin(1,1)+EA_hat/ds)/A_skin(1,1);
beta = sqrt(D_skin(2,2)/D_skin(1,1));
AR_bar = b/ds;
AR = b/a;
%% Section 9.2.1.1 + 9.2.1.2- 
Force = 7*10^3;
F_skin = A_skin(1,1)/(A_skin(1,1)+ EA_hat/ds)*Force; % eq 9.18
F_hat = Force - F_skin;
N_xskin = F_skin/a; % eq 9.19
m = 2; % makes realistic a_eff - first condition of of skin stiffener seperation
k = 10;
N_x_skin = (pi.^2./b.^2).*(D_skin(1,1).*k.^2 + 2.*(D_skin(2,1)+2.*D_skin(3,3)).*AR_bar.^2 + D_skin(2,2).*AR_bar.^4./k.^2); % eq 9.20
N_x_panel = (pi.^2./b.^2).*((D_skin(1,1)+EI_hat/ds).*m.^2 + 2.*(D_skin(2,1)+2.*(n_stiff*D_hat(3,3))).*AR.^2 + D_skin(2,2).*AR.^4./m.^2); % eq 9.20
%N_panel = (pi^2/a^2)*(D_eq(1,1)*m^2 + 2*(D_eq(2,1)+2*D_eq(3,3))*AR^2 + D_eq(2,2)*AR^4/m^2); % eq 9.20
a_eff = a/(2*(1+2*(1+A_skin(1,2)/A_skin(1,1))*(1-N_x_skin(1)/N_xskin)*(A_skin(1,1)/(A_skin(1,1)+3*A_skin(2,2)))));
a_eff = 0.303*a;
EI_skin_buckles = D_skin(1,1)*ds*((beta)*(2*lambda*AR_bar^2-beta*AR^4)+ 2*((D_skin(1,2)+2*D_skin(3,3))/D_skin(1,1))*(lambda*AR_bar^2-AR^2)-1); % Eq 9.30
k_opt = (D_skin(2,2)/(D_skin(1,1)))^0.25*(AR_bar);
m_opt = (D_skin(2,2)/(D_skin(1,1)+EI_hat/ds))^0.25*(AR); % m=1 since m<1 eq
k = 10;
N_x_skin = (pi.^2./b.^2).*(D_skin(1,1).*k.^2 + 2.*(D_skin(2,1)+2.*D_skin(3,3)).*AR_bar.^2 + D_skin(2,2).*AR_bar.^4./k.^2); % eq 9.20



skin_buckles = EI_hat > EI_skin_buckles;
disp(skin_buckles)
% if 1 then skin buckles before stiff


% stiffener buckling - Panel breaker condition - didnt need to do but ok 
F_skin_eff = 2*A_skin(1,1)*a_eff/(2*A_skin(1,1)*a_eff + EA_hat)*Force; %E 9.31
F_hat_ = EA_hat/(2*A_skin(1,1)*a_eff + EA_hat)*Force; %E 9.32
% Proper Deconstruction of distribution of force

F_hat_single = F_hat_*(ds/a); % Eq 9.33
F_hat_b = pi^2*EI_hat/b^2; %Eq 9.34
F_skin_bucking = a*N_x_skin;
k = k_opt;
% no post buckling PB = 1
F_skin_b = a*(pi^2/b^2)*(D_skin(1,1)*k^2 + 2*(D_skin(2,1)+2*D_skin(3,3))*AR_bar^2 + D_skin(2,2)*AR_bar^4/k^2);

%EI_hat_buckles = (lambda-1)*ds^2/(2*a_eff)*(D_skin(1,1)*k^2 + 2*(D_skin(2,1)+2*D_skin(3,3))*AR_bar^2 + D_eq(2,2)*AR_bar^4/k^2);
EI_hat_buckles = (lambda-1)*ds*(D_skin(1,1)*k^2 + 2*(D_skin(2,1)+2*D_skin(3,3))*AR_bar^2 + D_eq(2,2)*AR_bar^4/k^2);

stiff_buckles = EI_hat < EI_hat_buckles;
disp(stiff_buckles)
% Confirms that stiff do not buckles on skin does
% in the end the skin and stiffeners buckle with 

%% Crippling calculations
% hat
inv_A_skin = inv(A_skin);
E_skin = 1/(inv_A_skin(1,1)*h_skin);
m = 1;
N_0 = pi^2/b^2*(D_skin(1,1)*m^2 + 2*(D_skin(1,2) + 2*D_skin(3,3))*AR^2 + D_skin(2,2)*AR^4/m^2);
F_0 = N_0*a;

b_m1 = 1.5/100; % Flange 
b_m2 = 1.5/100; % Angle 
b_m3 = 2/100; % Top
ratio_m1 = b_m1/h_hat;
ratio_m2 = b_m2/h_hat;
ratio_m3 = b_m3/h_hat;
ratio_m4 = a_eff/h_skin;

% Stiffener
crip_m1 = OEFcrippling(ratio_m1,Y_C);
crip_m2 = NEFcrippling(ratio_m2,Y_C);
crip_m3 = NEFcrippling(ratio_m3,Y_C);
% Skin 
crip_m4 = NEFcrippling(ratio_m4,Y_C);


%%


%%
%F_skin = A_skin(1,1)/(A_skin(1,1) + A_hatT(1,1))*Force; %eq. (9.18)
%N_skin = F_skin/bp; %eq. (9.19)

% Assumed simply supported
%N_x_skin = (pi^2/a^2)*(D_skin(1,1)*k^2 + 2*(D_skin(1,2)+2*D_skin(3,3))*AR^2 + D_skin(2,2)*AR^4/k^2); %(9.20)

%m = 1;

%N_x_panel = (pi^2/a^2)*(D_eq(1,1)*m^2 + 2*(D_eq(2,1)+2*D_eq(3,3))*AR^2 + D_eq(2,2)*AR^4/m^2);





%% Functions
function [A,B,D,ABD,h_p,Qbar] = ABD_matrixCal(thetadt,E1,E2,nu12,G12)
Nplies_p  = length(thetadt);
thetadb = fliplr(thetadt); % ply angles in degrees, from bottom
h_ply   = 0.1905*10^(-3)/2;   % SI units, meters
h_p      = Nplies_p * h_ply ;
z       = -h_p/2:h_ply:h_p/2;

% Creation of layering position in ply starting from bottom
for i = 1:Nplies_p
    zbar(i) = - (h_p + h_ply)/2 + i*h_ply;
end

% ABD Matrix Initialisation
A  = zeros(3,3);
B  = zeros(3,3);
D  = zeros(3,3);

% Creation of Reduced Compliance and Stiffness Matrix
[S, Q] = ReducedComplianceStiffness(E1,E2,nu12,G12);
[moduli]= [E1 E2 nu12 G12];

% Creation of ABD Matrix
for i = 1:Nplies_p
    [Qbar(:,:,i),Sbar(:,:,i)] = QbarandSbar(thetadb(i),moduli);
    A  = A  + Qbar(:,:,i) * (z(i+1)-z(i)) ; %N/m,
    B  = B  + (1/2)* Qbar(:,:,i)* (z(i+1)^2-z(i)^2); %N
    D  = D  + (1/3)* Qbar(:,:,i) * (z(i+1)^3-z(i)^3); %Nm
    ABD  = [A  B ; B  D ];
end
end
% Take I beam stiffness with same thickness everywhere
% Recalculate ABD matrix




function [r] = NEFcrippling(ratio,Y_C)
r = Y_C*11/(ratio)^(1.124);
end
function [r] = OEFcrippling(ratio,Y_C)
r = Y_C*1.63/(ratio)^(0.717);
end

function [P_cr_n, x_axis] = BucklingElastic(EI,L,k,m)
% Assumed Pin Pin condition
P_cr_n = m.^2 + k.*L.^4./(pi.^4.*EI.*m^2);
x_axis = k.*L.^4./(pi.^4.*EI.*m^2);
end

function [P_cr] = Buckling(EI,L)
% Assumed Pin Pin condition
P_cr = pi.^2.*EI./L.^2;

end


function [FI_1,FI_2,FI_3]= MaxStress(sigma1,sigma2,sigma3,X_T,X_C,Y_T,Y_C,S_f)
% This function returns the Failure indices under the maximum stress criterion
% for fiber-reinforced materials.
% There are inputs are the local stresses with the materials strengths
% and returns each failure index
if sigma1>=0
    FI_1 = sigma1/X_T;
else
    FI_1 = -sigma1/X_C;
end

if sigma2>=0
    FI_2 = sigma2/Y_T;
else
    FI_2 = -sigma2/Y_C;
end
FI_3 = abs(sigma3)/S_f;
end


function [Qbar,Sbar] = QbarandSbar(angle,moduli)
% This function returns the transformed reduced stiffness
% matrix for fiber-reinforced materials.
% There are four arguments representing four
% material constants. The size of the reduced
% stiffness and compliance matrix is 3 x 3.
% Sine of the angle of the lamina
s=sind(angle);
% Cosine of the angle of the lamina
c=cosd(angle);
% transformed compliance matrix values 
Sb11=((1/moduli(1))*c^4)+((2*(-moduli(3)/moduli(1))+(1/moduli(4)))*s^2*c^2)+((1/moduli(2))*s^4);
Sb12=((-moduli(3)/moduli(1))*(s^4+c^4))+(((1/moduli(1))+(1/moduli(2))-(1/moduli(4)))*s^2*c^2);
Sb16=((2*(1/moduli(1))-2*(-moduli(3)/moduli(1))-(1/moduli(4)))*s*c^3)-((2*(1/moduli(2))-2*(-moduli(3)/moduli(1))-(1/moduli(4)))*s^3*c);
Sb22=((1/moduli(1))*s^4)+((2*(-moduli(3)/moduli(1))+(1/moduli(4)))*s^2*c^2)+((1/moduli(2))*c^4);
Sb26=((2*(1/moduli(1))-2*(-moduli(3)/moduli(1))-(1/moduli(4)))*s^3*c)-((2*(1/moduli(2))-2*(-moduli(3)/moduli(1))-(1/moduli(4)))*s*c^3);
Sb66=(2*(2*(1/moduli(1))+2*(1/moduli(2))-4*(-moduli(3)/moduli(1))-(1/moduli(4)))*s^2*c^2)+((1/moduli(4))*(s^4+c^4));
% transformed compliance matrix and transformed reduced stiffness matrix
Sbar=[Sb11 Sb12 Sb16;Sb12 Sb22 Sb26;Sb16 Sb26 Sb66];
Qbar=inv(Sbar);
end


function [S, Q] = ReducedComplianceStiffness(E_1,E_2,v_12,G_12)
% This function returns the reduced stiffness
% matrix for fiber-reinforced materials.
% There are four arguments representing four
% material constants. The size of the reduced
% stiffness and compliance matrix is 3 x 3.
v_21 = v_12*E_2/E_1;
S = [1/E_1 -v_12/E_1 0 ;
    -v_12/E_1 1/E_2 0;
    0 0 1/G_12];
Q = [E_1/(1-v_12*v_21), v_12*E_2/(1-v_12*v_21), 0 ;
    v_12*E_2/(1-v_12*v_21), E_2/(1-v_12*v_21), 0 ; 
    0,0,G_12];

end

function [strain_loc] = strain_gtol(strain_glo,angle)
% Sine of the angle of the lamina
s=sind(angle);
% Cosine of the angle of the lamina
c=cosd(angle);
% Reuter Matrix
R=[1 0 0; 0  1  0;0  0  2];
% Transformation matrix
T=[c^2, s^2, 2*s*c; s^2, c^2, -2*s*c; -s*c, s*c, c^2-s^2;];
% local strain
strain_loc=R*T*inv(R)*strain_glo;
%strain_loc=T*strain_glo;

end

function [stress_loc] = stress_gtol(stress_glo,angle)
% Global stresses
% Sine of the angle of the lamina
s=sind(angle);
% Cosine of the angle of the lamina
c=cosd(angle);
% Transformation matrix
R = [1 0 0; 0  1  0;0  0  2];
% Transformation matrix
T = [c^2, s^2, 2*s*c; s^2, c^2, -2*s*c; -s*c, s*c, c^2-s^2;];
% local strain
stress_loc=T*stress_glo;
end