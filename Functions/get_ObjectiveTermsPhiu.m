function [Hu,gu,rhou] = get_ObjectiveTermsPhiu(Qu, Nph, Uk_bar)
% Compute the objective terms for Phiu

Wu = chol(Qu,'lower');
I_N = eye(Nph);
Wu_bar = kron(I_N,Wu);
Hu = Wu_bar'*Wu_bar;
Mu = -Hu;
gu = Mu*Uk_bar;
rhou = (1/2)*Uk_bar'*Wu_bar'*Wu_bar*Uk_bar;

end