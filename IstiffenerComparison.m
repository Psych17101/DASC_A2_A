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
n_stiff=4; %number of stiffeners
ds = a/n_stiff; 
%% Initialisation of Composites Model
% ABD Skin/Panel
% Laminate definition (plies of equal thickness)
%thetadt_skin = [0 90 +45 -45 -45 +45 90 0];%Original 
thetadt_skin = [+45 -45 0 0 45 +45 ];% New tested
%thetadt_skin = [0 90 0 0 90 0];% New tested
[A_skin,B_skin,D_skin,ABD_skin,h_skin,Qbar] = ABD_matrixCal(thetadt_skin,E1,E2,nu12,G12);
% Calcuation of global strain for first ply failure

F = [-7e3/a;0;0];
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

% ABD I stiffeners
thetadt_I1 = [+45 -45 90 0 0 0 0 90-45 +45];% bottom
[A_I1,B_I1,D_I1,ABD_I1,h_I1,Qbar_1] = ABD_matrixCal(thetadt_I1,E1,E2,nu12,G12);
thetadt_I2 = [+45 -45 90 0 0 90 -45 +45];% web
[A_I2,B_I2,D_I2,ABD_I2,h_I2,Qbar_2] = ABD_matrixCal(thetadt_I2,E1,E2,nu12,G12);
thetadt_I3 = [+45 -45 90 0 0 90 -45 +45];% Top
[A_I3,B_I3,D_I3,ABD_I3,h_I3,Qbar_3] = ABD_matrixCal(thetadt_I3,E1,E2,nu12,G12);

if h_I1/h_skin >= 1.5 ||  h_I1/h_skin < 1
    fprintf("Thicknes of flange in norms!\n")
end

%% Properties of Stiffeners
%I
Area_I = 8.2078e-05;%Cross-sectional area of hat stiffener [m^2]
Interia_I =  2.0394e-08;
members_I = 3;
b_I1 = 1.5/100;
b_I2 = 2/100;
b_I3 = 1.5/100; % For t stiffener

%E_
inv_A_I1 = inv(A_I1);
inv_A_I2 = inv(A_I2);
inv_A_I3 = inv(A_I3);
E_I1 = 1/(inv_A_I1(1,1)*h_I1);
E_I2 = 1/(inv_A_I2(1,1)*h_I2);
E_I3 = 1/(inv_A_I3(1,1)*h_I3);
Area_I1 = h_I1*b_I1;
Area_I2 = h_I2*b_I2;
Area_I3 = h_I3*b_I3;

EA_I = E_I1*Area_I1 + E_I2*Area_I2 + E_I3*Area_I3;

inv_D_I1 = inv(D_I1);
inv_D_I2 = inv(D_I2);
inv_D_I3 = inv(D_I3);
Eb_I1 = 12/(inv_D_I1(1,1)*h_I1^3);
Eb_I2 = 12/(inv_D_I2(1,1)*h_I2^3);
Eb_I3 = 12/(inv_D_I3(1,1)*h_I3^3);
ACy = 0.011805;
d_I1 = ACy - h_I1/2;
d_I3 = b_I2 + h_I1 + h_I3/2 - ACy;
%d_I3 = 0 % for T stiffener
EI_I = Eb_I1*(b_I1*h_I1^3/12 + Area_I1*d_I1^2) + Eb_I2*h_I2*b_I3^3/12 + Eb_I3*(b_I3*h_I3^3/12 + Area_I3*d_I3^2);

lambda = (A_skin(1,1)+EA_I/ds)/A_skin(1,1);
beta = sqrt(D_skin(2,2)/D_skin(1,1));
AR_bar = b/ds;
AR = b/a;

%% Section 9.2.1.1 + 9.2.1.2
Force = 7*10^3;
F_skin = A_skin(1,1)/(A_skin(1,1)+ EA_I/ds)*Force; % eq 9.18
N_xskin = F_skin/a; % eq 9.19
m = 2; % makes realistic b_eff - first condition of of skin stiffener seperation
N_x_skin = (pi.^2./b.^2).*(D_skin(1,1).*m.^2 + 2.*(D_skin(2,1)+2.*D_skin(3,3)).*AR_bar.^2 + D_skin(2,2).*AR_bar.^4./m.^2); % eq 9.20
%N_panel = (pi^2/a^2)*(D_eq(1,1)*m^2 + 2*(D_eq(2,1)+2*D_eq(3,3))*AR^2 + D_eq(2,2)*AR^4/m^2); % eq 9.20
b_eff = a/(2*(1+2*(1+A_skin(1,2)/A_skin(1,1))*(1-N_x_skin/N_xskin)*(A_skin(1,1)/(A_skin(1,1)+3*A_skin(2,2)))));
b_eff = 0.303*a;
EI_skin_buckles = D_skin(1,1)*ds*((beta)*(2*lambda*AR_bar^2-beta*AR^4)+ 2*(D_skin(1,2)+2*D_skin(3,3))/D_skin(1,1)*(lambda*AR_bar^2-AR^2)-1); % Eq 9.30
k_opt = (D_skin(2,2)/(D_skin(1,1)))^0.25*(AR_bar);
m_opt = (D_skin(2,2)/(D_skin(1,1)+EI_I/ds))^0.25*(AR); % m=1 since m<1 eq

skin_buckles = EI_I > EI_skin_buckles;
%disp(skin_buckles)
% if 1 then skin buckles before stiff
skin_buckles_forces = N_xskin > N_x_skin;
if skin_buckles_forces == 1
    fprintf("Skin buckling happens!\n")
    fprintf("Force in skin = %d \n", N_xskin )
    fprintf("skin buckling force = %d \n", N_x_skin )
end
disp(skin_buckles_forces)

F_skin_eff = 2*A_skin(1,1)*b_eff/(2*A_skin(1,1)*b_eff + EA_I)*Force; %E 9.31
F_I_ = EA_I/(2*A_skin(1,1)*b_eff + EA_I)*Force; %E 9.32
F_I_1 = F_I_*(ds/a); % Eq 9.33
F_I_b = pi^2*EI_I/a^2; %Eq 9.34
F_max_I_b = F_I_b/F_I_1*F_I_;

stiff_buckles_forces = F_I_1 > F_I_b;
if stiff_buckles_forces == 1
    fprintf("Skin buckling happens!'n")
end
disp(stiff_buckles_forces)


k = k_opt;
% no post buckling PB = 1
F_skin_b = b*(pi^2/a^2)*(D_skin(1,1)*k^2 + 2*(D_skin(2,1)+2*D_skin(3,3))*AR_bar^2 + D_skin(2,2)*AR_bar^4/k^2);

%EI_I_buckles = (lambda-1)*ds^2/(2*b_eff)*(D_skin(1,1)*k^2 + 2*(D_skin(2,1)+2*D_skin(3,3))*AR_bar^2 + D_skin(2,2)*AR_bar^4/k^2);
EI_I_buckles = (lambda-1)*ds*(D_skin(1,1)*k^2 + 2*(D_skin(2,1)+2*D_skin(3,3))*AR_bar^2 + D_skin(2,2)*AR_bar^4/k^2);

stiff_buckles = EI_I < EI_I_buckles;
%disp(stiff_buckles)
% Confirms that stiff do not buckles on skin does
% in the end the skin and stiffeners buckle with 

%% Crippling calculations
% hat
ratio_m1 = b_I1/h_I1;
ratio_m2 = b_I2/h_I2;
ratio_m3 = b_I3/h_I3;
ratio_m4 = b_eff/(h_skin+h_I1);

% Stiffener
crip_m1 = OEFcrippling(ratio_m1,Y_C);
crip_m2 = NEFcrippling(ratio_m2,Y_C);
crip_m3 = NEFcrippling(ratio_m3,Y_C);
% Skin 
crip_m4 = NEFcrippling(ratio_m4,Y_C);

%%
inv_A_skin = inv(A_skin);
E_skin = 1/(inv_A_skin(1,1)*h_skin);
m = 1;
N_0 = pi^2/b^2*(D_skin(1,1)*m^2 + 2*(D_skin(1,2) + 2*D_skin(3,3))*AR^2 + D_skin(2,2)*AR^4/m^2);
F_0 = N_0*a;


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