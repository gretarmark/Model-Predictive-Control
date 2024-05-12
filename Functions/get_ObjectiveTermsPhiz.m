function [Hz,gz,rhoz] = get_ObjectiveTermsPhiz(Qe, Nph, Gamma,rk,bk)
% Compute the objective terms for Phiz

Qz = Qe;
Wz = chol(Qz,'lower');
I_N = eye(Nph);
Wz_bar = kron(I_N,Wz);
WzG = (Wz_bar*Gamma);
Mz = -WzG'*Wz_bar;
Hz = WzG'*WzG;  
ck = rk-bk;
gz = Mz*ck;
rhoz = (1/2)*ck'*Wz_bar'*Wz_bar*ck;

end