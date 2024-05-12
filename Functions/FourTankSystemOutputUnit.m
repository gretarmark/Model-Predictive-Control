function z = FourTankSystemOutputUnit(x,p)

rho = p(12);
A = p(5:8);
z = zeros(2,1);
z1 = x(1,1)/(rho*A(1,1));
z2 = x(2,1)/(rho*A(2,1));
z = [z1; z2]; 
end