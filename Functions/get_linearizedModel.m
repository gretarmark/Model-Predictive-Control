function [Ac,Bc,Cc,Cz,Ec] = get_linearizedModel(p,xs)

a1 = p(1,1);
a2 = p(2,1);
a3 = p(3,1);
a4 = p(4,1);

A1 = p(5,1);
A2 = p(6,1);
A3 = p(7,1);
A4 = p(8,1);

gamma1  = p(9,1);
gamma2  = p(10,1);
g       = p(11,1);
rho     = p(12,1);

syms t u_1 u_2 x_1(t) x_2(t) x_3(t) x_4(t) F3 F4 d_1 d_2

x10 = xs(1);
x20 = xs(2);
x30 = xs(3);
x40 = xs(4);

Dx_1 = diff(x_1);
Dx_2 = diff(x_2);
Dx_3 = diff(x_3);
Dx_4 = diff(x_4);

eq1 = rho*gamma1*u_1 + rho*a3*sqrt(2*g*(x_3/(rho*A3))) - rho*a1*sqrt(2*g*(x_1/(rho*A1))) == Dx_1;
eq2 = rho*gamma2*u_2 + rho*a4*sqrt(2*g*(x_4/(rho*A4))) - rho*a2*sqrt(2*g*(x_2/(rho*A2))) == Dx_2;
eq3 = rho*(1 - gamma2)*u_2 - rho*a3*sqrt(2*g*(x_3/(rho*A3))) + rho*d_1 == Dx_3;
eq4 = rho*(1 - gamma1)*u_1 - rho*a4*(sqrt(2*g*(x_4/(rho*A4)))) + rho*d_2 == Dx_4;
[V,Y] = odeToVectorField([eq1;eq2;eq3;eq4;]);

V = strrep(string(V),'Y[1]',string(Y(1)));
V = strrep(string(V),'Y[2]',string(Y(2)));
V = strrep(string(V),'Y[3]',string(Y(3)));
V = strrep(string(V),'Y[4]',string(Y(4)));
    
V = [V(2) ; V(1) ; V(3) ; V(4)];
V = str2sym(V);
        
A = (jacobian(V,Y));
A = vpa([A(:,2) A(:,1) A(:,3) A(:,4)],8);
A = subs(subs(subs(subs(A,x_1,xs(1)),x_2,xs(2)),x_3,xs(3)),x_4,xs(4));
% A = subs(subs(subs(subs(A,x_1,ys(1)),x_2,ys(2)),x_3,ys(3)),x_4,ys(4));
A = vpa(A,7);
Ac = cast(A,'double');
    
B = (jacobian(V,[u_1 ; u_2]));
B = vpa(B,8);
Bc = cast(B,'double');
    
At = [A1 ; A2 ; A3 ; A4];
Cc=diag(1./(rho*At));
Cz=Cc(1:2,:);

E = jacobian(V,[d_1 ; d_2]);
E = vpa(E,4);
Ec = cast(E,'double');

end