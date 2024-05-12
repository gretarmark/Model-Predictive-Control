function [p,a,A,g,rho] = get_parameters(gamma1,gamma2)

% Parameters
a1 = 1.2272;     %[cm2] Area of outlet pipe 1
a2 = 1.2272;     %[cm2] Area of outlet pipe 2
a3 = 1.2272;     %[cm2] Area of outlet pipe 3
a4 = 1.2272;     %[cm2] Area of outlet pipe 4

A1 = 380.1327;   %[cm2] Cross sectional area of tank 1
A2 = 380.1327;   %[cm2] Cross sectional area of tank 2
A3 = 380.1327;   %[cm2] Cross sectional area of tank 3
A4 = 380.1327;   %[cm2] Cross sectional area of tank 4

g = 981;        %[cm/s2] The acceleration of gravity
rho = 1.00;     %[g/cm3] Density of water

p = [a1; a2; a3; a4; A1; A2; A3; A4; gamma1; gamma2; g; rho];

a = p(1:4,1)';
A = p(5:8,1)';
end