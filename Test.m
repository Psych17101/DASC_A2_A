
a = 480e-3;
b_eff = 0.303*a;
AR = 480e-3;
D_11 = 659.7e-3;
D_12 = 466.9e-3;
D_22 = 659.7e-3;
D_33 = 494.0e-3;
m = 2;
N_x_skin = (pi^2/b_eff^2)*(D_11*m^2 + 2*(D_12+2*D_33)*AR^2 + D_22*AR^4/m^2); %(6.7)

%N/m