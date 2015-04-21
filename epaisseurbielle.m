function [ output_args ] = untitled( input_args )
MaxPressure= 6.76 ;
CompressionRatio = 8.2;
LOR=128.85;
PDiam=75;
Stroke=73.4;
FactorOfSav=1;
a=1/7500;
W=pi*(PDiam)^2*MaxPressure;
Wb=W*FactorOfSav
A = 11*410;
B = -Wb;
C = -Wb * (LOR/1.78) * a;
rho = B*B - 4*A*C;
t = sqrt((-B + sqrt(rho)) / (2*A))
end

