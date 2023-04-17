clear; format;
close;
format short g;

% Define material properties
Ex = 62.046*10^9;
Ey = 62.046*10^9;
vxy = 0.05 ;
Gxy  = 6*10^9  ; % Pa
X_T = 1517*10^6;
X_C = 1379*10^6;
Y_T = 1450*10^6;
Y_C = 1379*10^6;
S = 99*10^6;

% Panel dimensions
a = 1; % [m] y direction 
b = 0.48; % [m] x direction (side where force is applied)
n_S=5; %number of stiffeners
ds = b/n_S; 
P=7e3;
n_IS=6;
ds_new=b/n_IS;
%% base values

%layups of skin and hat and new skin
ogskin=[0 90 +45 -45 -45 +45 90 0];
layupskin=[0 90 +45 -45 -45 +45 90 0];
%layupskin=[ +45 -45 0 0 0 0 -45 +45];
oghat = [0 90 0 90 90 0 90 0];
hat=[0 90 0 0 90 0];
%hat=[+45 -45 0 0 0 0 -45 +45];

layupskin_45 = [ +45 -45 0 0 0 0 -45 +45];
%layupskin_45 = [ 0 90 0 0 0 0 90 0];
%layupskin_45=[0 90 +45 -45 -45 +45 90 0];

[E_og_skin, EA_og_skin, area_og_skin,A_og_skin,Qbar_og_skin]=calculateEAskin(ogskin, Ex, Ey, vxy, Gxy, b);
[E_skin, EA_skin, area_skin,A_skin,Qbar_skin]=calculateEAskin(layupskin, Ex, Ey, vxy, Gxy, b);
[E_skin_45, EA_skin_45, area_skin_45,A_skin_45,Qbar_skin_45]=calculateEAskin(layupskin_45, Ex, Ey, vxy, Gxy, b);


[E_og_stiff, EA_og_stiff, EI_og_stiff, area_og_stiff,A_og_hat,Qbar_og_hat]=calculateEAhat(oghat, Ex, Ey, vxy, Gxy, n_S, 1);
[E_stiff,EA_stiff, EI_stiff,area_stiff,A_hat,Qbar_hat]=calculateEAhat(hat, Ex, Ey, vxy, Gxy, n_IS, 0.8);

[E_I,EA_I, EI_I,area_I] = IStiff_prop(n_IS);

%% localSkinBuckling, input(n_S, layupskin), output(PB)

lr_og=EA_og_skin/(EA_og_skin+EA_og_stiff)
lr_new=EA_skin/(EA_skin+EA_stiff)
lr_I=EA_skin_45/(EA_skin_45+EA_I)

PB_og=localSkinBuckling(n_S, ogskin, EA_og_skin, EA_og_stiff, P)
PB_new=localSkinBuckling(n_IS, layupskin, EA_skin, EA_stiff, P)
PB_I=localSkinBuckling(n_IS, layupskin_45, EA_skin, EA_I, P)

Pcr=stiffenerColumnBuckling(EI_og_stiff);
P_og_skin = P*(1-lr_og);

CB_og = P_og_skin/Pcr

Pcr=stiffenerColumnBuckling(EI_stiff);
P_skin = P*(1-lr_new);

CB_new = P_skin/Pcr

Pcr=stiffenerColumnBuckling(EI_I);
P_skin = P*(1-lr_I);

CB_I = P_skin/Pcr


wr_hat=((area_og_skin+area_og_stiff)/(area_skin+area_stiff))^-1
wr_I=((area_og_skin+area_og_stiff)/(area_skin_45+area_I))^-1

%% Failure 
FI_max_og_skin = Failure_Criterion(P_og_skin,ogskin,A_og_skin,Qbar_og_skin(:,:,1),X_T,X_C,Y_T,Y_C,S)
FI_max_skin = Failure_Criterion(P_skin,layupskin,A_skin,Qbar_skin(:,:,1),X_T,X_C,Y_T,Y_C,S)
FI_max_og_hat = Failure_Criterion(P-P_og_skin,oghat,A_og_hat,Qbar_og_hat(:,:,1),X_T,X_C,Y_T,Y_C,S)
FI_max_hat = Failure_Criterion(P-P_skin,hat,A_hat,Qbar_hat(:,:,1),X_T,X_C,Y_T,Y_C,S)

%% plot
% 
% f=0.4:0.001:1;
% y1=[];
% y2=[];
% y3=[];
% y4=[];
% for i=f
%     
%     [E_stiff,EA_stiff, EI_stiff,area_stiff]=calculateEAhat(hat, Ex, Ey, vxy, Gxy, n_IS, i);
%     y1=[y1, EA_skin/(EA_skin+EA_stiff)]; %lr
%     y2=[y2, localSkinBuckling(n_IS, layupskin, EA_skin, EA_stiff, P)]; %PB
%     y3=[y3, P*(1- EA_skin/(EA_skin+EA_stiff))/stiffenerColumnBuckling(EI_stiff)];%CB
%     y4=[y4, ((area_og_skin+area_og_stiff)/(area_skin+area_stiff))^-1];
% 
% end
% plot(f, y1, f,y2,f,y3,f,y4, f, ones(size(f)), 0.708, 0.9191, "*", 0.8, y2(401), "*", 0.652, y3(253), "*"); legend("loading ratio", "Local buckling condition", "Column buckling condition", "Weight ratio", "Failure line"); xlabel("factor"); grid on;
% 
% [E_stiff,EA_stiff, EI_stiff,area_stiff]=IStiff_prop(n_IS); %there are now 6 stiffeners
% [E_stiff,EA_stiff, EI_stiff,area_stiff]=SStiffener_comp(n_IS);
% [E_stiff,EA_stiff, EI_stiff,area_stiff]=TStiffener_comp(n_IS);  
% 
% 
% %401
% %253
% 
% 
% % from the plot, the maximum 











%% Functions
function [Pcr]=stiffenerColumnBuckling(EI)
%EI=115.42;
L=1;
Pcr=pi^2*EI/L;
end
%failure!


function [E, EA, Area,A,Qbar] = calculateEAskin(layup, Ex, Ey, vxy, Gxy, b)
thetadt=layup;
[A,B,D,ABD,h_p,Qbar] = ABD_matrixCal(thetadt,Ex,Ey,vxy,Gxy);
aa=inv(A);
E10=1/aa(1,1)/h_p;
E20=1/aa(2,2)/h_p;
G12=1/aa(3,3)/h_p;
EA=E10*b*h_p;
Area=b*h_p;
E=E10;
end


function [E, EA, EI, Area,A,Qbar] = calculateEAhat(layup, Ex, Ey, vxy, Gxy, n_S, factor1) %for the original stiffener
thetadt=layup;
[A,B,D,ABD,h_p,Qbar] = ABD_matrixCal(thetadt,Ex,Ey,vxy,Gxy);
aa=inv(A);
E10=1/aa(1,1)/h_p;
E20=1/aa(2,2)/h_p;
G12=1/aa(3,3)/h_p;
w_stiff=1e-3.*[15, 15, 20, 15, 15]';
w_stiff(2)=factor1.*w_stiff(2);
w_stiff(4)=factor1.*w_stiff(4);
w_stiff(1)=factor1.*w_stiff(1);
w_stiff(end)=factor1.*w_stiff(end);
w_stiff(3)=factor1.*w_stiff(3);
w_stiff;
t_stiff=h_p.*ones(5,1);
E_stiff=E10.*ones(5,1);
Area_stiff=w_stiff.*t_stiff;
EA=n_S*sum(E_stiff.*Area_stiff);
Area=n_S*sum(Area_stiff);
E=E10;

phi=asind(0.5/(100*w_stiff(2)));
d_stiff=inv(D);
E_b =12./t_stiff.^3/d_stiff(1,1);
y=1e-2.*[t_stiff(1)/2, 1/2*1.5*cosd(phi), 1.5*cosd(phi)+t_stiff(3)/2, 1/2*1.5*cosd(phi), t_stiff(1)/2]';
y_neutral=sum(E_stiff.*Area_stiff.*y)/sum(E_stiff.*Area_stiff);
d=y_neutral-y;
EIind=E_b.*(w_stiff.*t_stiff.^3/12+Area_stiff.*d.^2);
EI=n_S*sum(EIind);


end


function [A,B,D,ABD,h_p,Qbar] = ABD_matrixCal(layup,E1,E2,nu12,G12)
thetadt=layup;
Nplies_p  = length(thetadt);
thetadb = fliplr(thetadt); % ply angles in degrees, from bottom
h_ply   = 0.1905*10^(-3)/2;   % SI units, meters
h_p      = Nplies_p * h_ply;
z       = -h_p/2:h_ply:h_p/2;

% ABD Matrix Initialisation
A  = zeros(3,3);
B  = zeros(3,3);
D  = zeros(3,3);

[moduli]= [E1 E2 nu12 G12];

% Creation of ABD Matrix
for i = 1:Nplies_p
    [Qbar(:,:,i),Sbar(:,:,i)] = QbarandSbar(thetadt(i),moduli);
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

% function [FI_max]=fpf(layup, F)
% % ABD Skin/Panel
% % Laminate definition (plies of equal thickness)
% %thetadt_skin = [0 90 +45 -45 -45 +45 90 0];%Original 
% thetadt_skin =layup;% New tested
% [A_skin,B_skin,D_skin,ABD_skin,h_skin,Qbar] = ABD_matrixCal(thetadt_skin,Ex,Ey,nuxy,Gxy);
% % Calcuation of global strain for first ply failure
% 
% F = [-F/a;0;0];
% strain_glo = A_skin\F;
% % Calculation of global stresses
% stress_glo = Qbar(:,:,1)*strain_glo; % global sigmaxx etc...
% Nplies = length(thetadt_skin);
% for l = 1:Nplies
%             % Calculations of local Strains and Stresses
%             [eps_loc] = strain_gtol(strain_glo,thetadt_skin(l));
%             [sigma_loc] = stress_gtol(stress_glo,thetadt_skin(l));% ply i angle in radians, from bottom
%             % Failure index with Maximum Stress criterion
%             [FI_1(l),FI_2(l),FI_3(l)]= MaxStress(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,S);
% end
% end