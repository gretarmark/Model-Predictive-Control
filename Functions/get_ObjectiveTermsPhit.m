function [Ht,gt] = get_ObjectiveTermsPhit(Qt2, Qt1, Nph, B)
% Compute the objective terms for Phit
% The objective funtion term related to the upper bound constraints

nu = size(B,2); % Inputs to the system

I_N = eye(Nph);
W_t2 = chol(Qt2,'lower');
W_t1 = chol(Qt1,'lower');
Wbar_t2 = kron(I_N,W_t2);
Wbar_t1 = kron(I_N,W_t1);
Ht = Wbar_t2'*Wbar_t2;
e = ones(1,nu*Nph)';
gt = Wbar_t1*e;

end