function [Hs,gs] = get_ObjectiveTermsPhis(Qs2,Qs1, Nph, B)
% Compute the objective terms for Phis
% The objective funtion term related to the lower bound constraints

nu = size(B,2); % Inputs to the system

I_N = eye(Nph);
W_s2 = chol(Qs2,'lower');
W_s1 = chol(Qs1,'lower');
Wbar_s2 = kron(I_N,W_s2);
Wbar_s1 = kron(I_N,W_s1);
Hs = Wbar_s2'*Wbar_s2;
e = ones(1,nu*Nph)';
gs = Wbar_s1*e;

end