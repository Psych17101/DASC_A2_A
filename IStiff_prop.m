function [E_I,EA_I, EI_I,Area_tot] = IStiffener_comp(n_s)
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

% ABD I stiffeners
thetadt_I1 = [0 90 0 0 90 0];% bottom
[A_I1,B_I1,D_I1,ABD_I1,h_I1,Qbar_1] = ABD_matrixCal(thetadt_I1,E1,E2,nu12,G12);
thetadt_I2 = [0 90 +45 -45 -45 +45 90 0];% web
[A_I2,B_I2,D_I2,ABD_I2,h_I2,Qbar_2] = ABD_matrixCal(thetadt_I2,E1,E2,nu12,G12);
thetadt_I3 = [0 90 0 0 0 0 0 0 90 0];% Top
[A_I3,B_I3,D_I3,ABD_I3,h_I3,Qbar_3] = ABD_matrixCal(thetadt_I3,E1,E2,nu12,G12);

b_I1 = 1/100;
b_I2 = 3/100;
b_I3 = 1/100; % For I stiffener

%[E_stiffI,EA_stiffI,EI_stiffI,Area_totI] = IStiffener_comp(A_I1,A_I2,A_I3,D_I1,D_I2,D_I3,h_I1,h_I2,h_I3,b_I1,b_I2,b_I3,n_s)


inv_A_I1 = inv(A_I1);
inv_A_I2 = inv(A_I2);
inv_A_I3 = inv(A_I3);
E_I1 = 1/(inv_A_I1(1,1)*h_I1);
E_I2 = 1/(inv_A_I2(1,1)*h_I2);
E_I3 = 1/(inv_A_I3(1,1)*h_I3);
E_I = E_I1+ E_I2 + E_I3;
Area_I1 = h_I1*b_I1;
Area_I2 = h_I2*b_I2;
Area_I3 = h_I3*b_I3;
EA_I1 = E_I1*Area_I1;
EA_I2 = E_I2*Area_I2;
EA_I3 = E_I3*Area_I3;
y1 = h_I1/2;
y2 = b_I2/2 + h_I1;
y3 = h_I1 + b_I2 + h_I3/2;
sumEAy = EA_I1*y1 + EA_I2*y2 + EA_I3*y3;
EA_I = EA_I1 + EA_I2 + EA_I3;

y_neutral = sumEAy/EA_I; 

inv_D_I1 = inv(D_I1);
inv_D_I2 = inv(D_I2);
inv_D_I3 = inv(D_I3);
Eb_I1 = 12/(inv_D_I1(1,1)*h_I1^3);
Eb_I2 = 12/(inv_D_I2(1,1)*h_I2^3);
Eb_I3 = 12/(inv_D_I3(1,1)*h_I3^3);
d_I1 = y_neutral - b_I2/2 - h_I1/2;
d_I3 = y_neutral + b_I2/2 + h_I1/2;
EI_I = Eb_I1*(b_I1*h_I1^3/12 + Area_I1*d_I1^2) + Eb_I2*h_I2*b_I3^3/12 + Eb_I3*(b_I3*h_I3^3/12 + Area_I3*d_I3^2);

Area_single_tot = Area_I1 + Area_I2 + Area_I3;
Area_tot = n_s*Area_single_tot;

end

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
