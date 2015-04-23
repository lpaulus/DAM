function bielle
R_e = 410;
E = 170000;
S = 2.2;
F_b = 4.8*10^4;
l_f = 0.1285;
lambda_l = pi*sqrt(E/R_e);
sigma = R_e/S;
I = 5.8563*4.94^4;

B = -F_b/sigma;
C = (F_b * (l_f)^2)/(sigma*lambda_l *I)-1;
A = B/C 

end