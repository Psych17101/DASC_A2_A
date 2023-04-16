function [FI_max] = Failure_Criterion(F,layup,A,Qbar,X_T,X_C,Y_T,Y_C,S)

F = [-F;0; 0];
strain_glo = A\F;
% Calculation of global stresses
stress_glo = Qbar*strain_glo; % global sigmaxx etc...

Nplies = length(layup);
% Calculations of Failure Criterion for each ply
for l = 1:Nplies
    [sigma_loc] = stress_gtol(stress_glo,layup(l));% ply i angle in radians, from bottom
    % Failure index with Maximum Stress criterion

    [FI_1(l),FI_2(l),FI_3(l)]= MaxStress(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,S);
end

FI_1_max = max(FI_1(:));
FI_2_max = max(FI_2(:));
FI_3_max = max(FI_3(:));
FI_max = max([FI_1_max,FI_2_max,FI_3_max]);

%% Functions
function [FI_1,FI_2,FI_3]= MaxStress(sigma1,sigma2,sigma3,X_T,X_C,Y_T,Y_C,S_f)

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
end

end
