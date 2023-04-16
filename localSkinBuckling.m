function [PB]=localSkinBuckling(n_S, layupskin, EA_skin, EA_stiff, P)

Ex=62.046e9;
Ey=Ex;
Gxy=4.826e9;
vxy=0.05;
vyx=vxy*Ey/Ex;
% theta=45;
% m=cosd(theta);
% n=sind(theta);
% 
% 
% 
% E1=(m^4/Ex+(1/Gxy-2*vxy/Ex)*m^2*n^2+n^4/Ey)^(-1);
% E2=(n^4/Ey+(1/Gxy-2*vxy/Ex)*m^2*n^2+m^4/Ex)^(-1);
% v12=E1*(vxy/Ex*(m^4+n^4)-(1/Ex+1/Ey-1/Gxy)*m^2*n^2);
% G12=(2*(2/Ex+2/Ey+4*vxy/Ex-1/Gxy)*m^2*n^2+1/Gxy*(m^2*n^2))^(-1);
lr=EA_skin/(EA_skin+EA_stiff);

% Q
Q=1-vxy*vyx;    
q=Q^-1;
Qxx=Ex*q; Qyy=Ey*q; Qxy=vxy*Ey*q; Qss=Gxy;
thetadt=layupskin;
[A,B,D,ABD,h_p,Qbar] = ABD_matrixCal(thetadt,Ex,Ey,vxy,Gxy);
aa=inv(A);
E10=1/aa(1,1)/h_p;
E20=1/aa(2,2)/h_p;
G12=1/aa(3,3)/h_p;

%m = 0:1:15; 
a = 1;
b=0.48;

ds=b/n_S;
mpred=a/ds*(D(2,2)/D(1,1))^(1/4); %-->12, 13
m=[floor(mpred), floor(mpred)+1];
%m=10;
AR_bar=a/ds;
N_x_skin = (pi^2/a^2)*(D(1,1).*m.^2 + 2*(D(1,2)+2*D(3,3))*AR_bar^2 + D(2,2)*AR_bar^4./m.^2); %(9.20) 
%plot(m, N_x_skin)
Force=lr*P;

PB=max(Force/b./N_x_skin); %1.13-->buckling, no good
end
%% 
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
