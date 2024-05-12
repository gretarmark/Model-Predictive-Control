function y = FourTankSystemSensor(x,p)

rho = p(12);
A = p(5:8);
y = zeros(4,1);
for i=1:4   
    y(i,1) = x(i,1)./(rho*A(i,1));
end
end