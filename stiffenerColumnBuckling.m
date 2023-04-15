function [Pcr]=stiffenerColumnBuckling(EI)
I=2.079e-9;
E=4.504e+10;
%EI=115.42;
L=1;
Pcr=pi^2*EI/L;
end
%failure!