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

a = 0.48; %m
l = 1.5/100; %m width of flange  


%% Initialisation of Composites Model
% ABD Stiffeners
thetadt_stiff = [0 90 0 90 90 0 90 0];% ply angles in degrees
[A_S, B_S,D_S, ABD_S, h_s] = ABD_matrixCal(thetadt_stiff,E1,E2,nu12,G12);

%% Calculations of membrane stiffness from eq 8.7 
% Current design - taking the approximations of b and thickness for each
% member of the hat stiffener
inv_A_S = inv(A_S);
E = 1/(inv_A_S(1,1)*h_s);
EA_stiff = 4*E*1.5/100*h_s + E*2/100*h_s;
fprintf('EA_stiff= ')
disp(EA_stiff)

inv_D_S = inv(D_S);
Eb = 12/(inv_D_S(1,1)*h_s^3);
A = h_s*l;
EI_stiff = 2*Eb*A*h_s^2/12 + (Eb*2/100*h_s^3/12 + A*0.009127^2) + 2*(Eb*A*h_s^2/12 + A*0.009127^2);
fprintf('EI_stiff= ')
disp(EI_stiff)

%Buckling Analysis
% normal column Buckling
P_cr_stiff = Buckling(EI_stiff, a);
fprintf('P_cr_stiff = ')
disp(P_cr_stiff)
% Buckling with elastic foundation (with panel)
k_range = linspace(0,X_T,11);
[P_cr_n_1, x_axis1] = BucklingElastic(EI_stiff,a,k_range,1);
[P_cr_n_2, x_axis2] = BucklingElastic(EI_stiff,a,k_range,2);
[P_cr_n_3, x_axis3] = BucklingElastic(EI_stiff,a,k_range,3);
%figure;
%hold on
%grid on
%plot(x_axis1,P_cr_n_1)
%plot(x_axis2,P_cr_n_2)
%plot(x_axis3,P_cr_n_3)

%Crippling Analysis
D_C = 2*D_S(3,3)+ 2*A_S(3,3)*h_s;
N_xcrit_flange = 12*D_C/a^2; % equ 6.13a % N/m

% estimation that 4kn applied on the stiffeners (distributed between skin and panel)
% so 4kn force over area of flange (1.5cm) = 2.6e-5 N/m

ratio_Lt = l/h_s; %Fig 8.18
opt_l = h_s*6;

ratio_crip1 = 2.151/(ratio_Lt)^0.717; % Equ 8.42
% B Basis
ratio_crip2 = 1.63/(ratio_Lt)^0.717; % Equ 8.42








%% changing layup
thetadt_stiff_new = [-45 45 90 0 0 90 45 -45];% ply angles in degrees, from top??? Not bottom?
[A_S_new , B_S_new, D_S_new, ABD_S_new, h_s_new] = ABD_matrixCal(thetadt_stiff_new,E1,E2,nu12,G12);

%% Calculations of membrane stiffness from eq 8.7
% Current design - taking the approximations of b and thickness for each
% member of the hat stiffener
inv_A_S_new = inv(A_S_new);
E_new = 1/(inv_A_S_new(1,1)*h_s);
EA_stiff_new = 4*E_new*1.5/100*h_s + E_new*2/100*h_s;
fprintf('EA_stiff_new = ')
disp(EA_stiff_new)

% With new lay - has higher membrane stiffness
inv_D_S_new = inv(D_S_new);
Eb_new = 12/(inv_D_S_new(1,1)*h_s^3);
EI_stiff_new = 2*Eb_new*A*h_s^2/12 + (Eb_new*2/100*h_s^3/12 + A*0.009127^2) + 2*(Eb_new*A*h_s^2/12 + A*0.009127^2);
fprintf('EI_stiff_new  = ')
disp(EI_stiff_new)

% Buckling Analysis
P_cr_stiff_new = Buckling(EI_stiff_new, a);
fprintf('P_cr_stiff = ')
disp(P_cr_stiff_new)

% Crippling Analysis
D_C_new = 2*D_S_new(3,3)+ 2*A_S_new(3,3)*h_s;
N_xcrit_flange_new = 12*D_C_new/a^2; % equ 6.13a
% Increased 


%% Calculations of membrane stiffness of figure 8.6 - but I beam (higher moment of inertia)
% Laminate definition (plies of equal thickness)
% Member 1
thetadt_I1  = [45 -45 90 0 0 90 90 0 0 90 -45 45];% ply angles in degrees, from top??? Not bottom?
[A_I1,B_I1,D_I1,ABD_I1,h_I1] = ABD_matrixCal(thetadt_I1,E1,E2,nu12,G12);

% Member 2
thetadt_I2  = [45 -45 -45 45 45 -45 -45 45];% ply angles in degrees, from top??? Not bottom?
[A_I2,B_I2,D_I2,ABD_I2,h_I2] = ABD_matrixCal(thetadt_I2,E1,E2,nu12,G12);

% Membre 3
thetadt_I3  = [45 -45 0 90 0 0 0 90 0-45 45];% ply angles in degrees, from top??? Not bottom?
[A_I3,B_I3,D_I3,ABD_I3,h_I3] = ABD_matrixCal(thetadt_I3,E1,E2,nu12,G12);

inv_A_I1 = inv(A_I1);
inv_A_I2 = inv(A_I2);
inv_A_I3 = inv(A_I3);

E_I1 = 1/(inv_A_I1(1,1)*h_I1);
E_I2 = 1/(inv_A_I2(1,1)*h_I2);
E_I3 = 1/(inv_A_I3(1,1)*h_I3);
EA_I = E_I1*h_I1*2/100 + E_I2*h_I2*3/100 + E_I3*h_I3*2/100;
fprintf('EA_I = ')
disp(EA_I)

inv_D_I1 = inv(D_I1);
inv_D_I2 = inv(D_I2);
inv_D_I3 = inv(D_I3);
Eb_I1 = 12/(inv_D_I1(1,1)*h_I1^3);
Eb_I2 = 12/(inv_D_I2(1,1)*h_I2^3);
Eb_I3 = 12/(inv_D_I3(1,1)*h_I3^3);
EI_I = Eb_I1*(2/100*h_I1^3/12 + 2/100*h_I1*(1.5/100)^2) + Eb_I2*(3/100*h_I3^2/12) + Eb_I3*(2/100*h_I2^3/12 + 2/100*h_I2*(1.5/100)^2);
fprintf('EI_I = ')
disp(EI_I)

% Buckling Analysis
P_cr_I = Buckling(EI_I, a);
fprintf('Pr_cr_I =')
disp(P_cr_I)


%% Functions
function [P_cr_n, x_axis] = BucklingElastic(EI,L,k,m)
% Assumed Pin Pin condition
P_cr_n = m.^2 + k.*L.^4./(pi.^4.*EI.*m^2);
x_axis = k.*L.^4./(pi.^4.*EI.*m^2);
end

function [P_cr] = Buckling(EI,L)
% Assumed Pin Pin condition
P_cr = pi.^2.*EI./L.^2;

end

function [A,B,D,ABD,h_p] = ABD_matrixCal(thetadt,E1,E2,nu12,G12)
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




