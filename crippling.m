
Ex=62.046e9;
Ey=Ex;
Gxy=4.826e9;
vxy=0.05;
vyx=vxy*Ey/Ex;
X_T = 1517*10^6;
X_C = 1379*10^6;
Y_T = 1450*10^6;
Y_C = 1379*10^6;
S = 99*10^6;
Force=7e3;
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


%% For skin
Q=1-vxy*vyx;    
q=Q^-1;
Qxx=Ex*q; Qyy=Ey*q; Qxy=vxy*Ey*q; Qss=Gxy;
thetadt=[0 90 45 -45 -45 45 90 0];
[A,B,D,ABD,h_p,Qbar] = ABD_matrixCal(thetadt,Ex,Ey,vxy,Gxy);
aa=inv(A);
E10=1/aa(1,1)/h_p;
E20=1/aa(2,2)/h_p;
G12=1/aa(3,3)/h_p;
a = 1;
b=0.48;
n_S=5;
ds=b/n_S;

%PBi=1/1.1364; %inverse of the PB=Px/Pcr, in our case

%b_eff = ds/(2*(1+2*(1+A(1,2)/A(1,1))*(1-PBi)*(A(1,1)/(A(1,1)+3*A(2,2)))));



%% For stiffener
Q=1-vxy*vyx;    
q=Q^-1;
Qxx=Ex*q; Qyy=Ey*q; Qxy=vxy*Ey*q; Qss=Gxy;
thetadt=[0 90 0 90 90 0 90 0];
[As,Bs,Ds,ABDs,h_ps,Qbars] = ABD_matrixCal(thetadt,Ex,Ey,vxy,Gxy);
aas=inv(As);
E10s=1/aas(1,1)/h_ps;
E20s=1/aas(2,2)/h_ps;
G12s=1/aas(3,3)/h_ps;



%% members
w_stiff=0.5*1e-3.*[15, 15, 20, 15, 15]';
t_stiff=h_ps.*ones(5,1);
E_stiff=E10s.*ones(5,1);
Area_stiff=w_stiff.*t_stiff;
EA_stiff=E_stiff.*Area_stiff;

% adding skin, do not add for pre buckling
w_skin=[0.48, 0.48]';
t_skin=[h_p, h_p]';
E_skin=[E10, E10]';
Area_skin=t_skin.*w_skin;
EA_skin=E_skin.*Area_skin;
loaddist=EA_skin(1)/(EA_skin(1)+6*sum(EA_stiff));
loaddistapprox=A(1,1)/(A(1,1)+sum(EA_stiff)/ds);

% total
%w=[w_stiff; w_skin];
%t=[t_stiff; t_skin];
%E=[E_stiff; E_skin];
%Area=[Area_stiff; Area_skin];
%EA=[EA_stiff; EA_skin];
Fishare=EA_stiff./sum(EA_stiff);

% members 1, 5 OEF rest NEF

func=@(x) x.^2/X_T/X_C+(1/X_T-1/X_C)-1;
SigmaCU=abs(fzero(func, 1));


Force_stiff=(1-loaddist)*Force;
AppliedF=Force_stiff.*Fishare;
AppliedSigma=AppliedF./Area_stiff;
ratiobt=w_stiff./t_stiff;
ratiocrtocu=[OEFcrippling(ratiobt(1)), NEFcrippling(ratiobt(2)), NEFcrippling(ratiobt(3)), NEFcrippling(ratiobt(4)), OEFcrippling(ratiobt(5))]';
failcrippling=SigmaCU.*ratiocrtocu;
allowedratio=AppliedSigma./failcrippling; % OK! it does not fail in crippling


% here we can see that the margin is quite high for the stiffener which
% means that we can transfer load from the skin, which fails in local
% buckling, to the stiffeners as they do not fail in cripppling








%% Stiffener Buckling, 
phi=19.471;
d_stiff=inv(Ds);
E_b =12./t_stiff.^3/d_stiff(1,1);
y=1e-2.*[t_stiff(1)/2, 1/2*1.5*cosd(phi), 1.5*cosd(phi)+t_stiff(3)/2, 1/2*1.5*cosd(phi), t_stiff(1)/2]';
y_neutral=sum(E_stiff.*Area_stiff.*y)/sum(E_stiff.*Area_stiff);
d=y_neutral-y;
EI=E_b.*(w_stiff.*t_stiff.^3/12+Area_stiff.*d.^2);
EItot=sum(EI);








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

function [r] = NEFcrippling(ratio)
r = 11/(ratio)^(1.124);
end
function [r] = OEFcrippling(ratio)
r = 1.63/(ratio)^(0.717);
end
