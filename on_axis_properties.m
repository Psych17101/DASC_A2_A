Ex=62.046;
Ey=Ex;
Gxy=4.826;
vxy=0.05;
theta=-45;
m=cosd(theta);
n=sind(theta);



E1=(m^4/Ex+(1/Gxy-2*vxy/Ex)*m^2*n^2+n^4/Ey)^(-1);
E2=(n^4/Ey+(1/Gxy-2*vxy/Ex)*m^2*n^2+m^4/Ex)^(-1);
v12=E1*(vxy/Ex*(m^4+n^4)-(1/Ex+1/Ey-1/Gxy)*m^2*n^2);
G12=(2*(2/Ex+2/Ey+4*vxy/Ex-1/Gxy)*m^2*n^2+1/Gxy*(m^2*n^2))^(-1);