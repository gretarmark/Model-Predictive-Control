function z = FourTankSystemOutput(x,p)

rho = p(12);
A = p(5:8);
z = zeros(2,1);
for i=1:2
    z(i,1) = x(i,1)./(rho*A(i,1));
end
end