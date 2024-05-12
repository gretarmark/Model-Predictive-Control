function y = FourTankSystemSensorUnit(x,p)

rho = p(12);
A = p(5:8);
y = zeros(4,1);
y(:,1) = x(:,1)./(rho*A(:,1));

end