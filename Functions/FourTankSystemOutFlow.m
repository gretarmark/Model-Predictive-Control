function Qoutflow = FourTankSystemOutFlow(h,p)
% y = h = height
a = p(1:4);
g = p(11);
Qoutflow = zeros(4,1);
Qoutflow(:,1) = a.*sqrt(2*g*h(:,1));

end