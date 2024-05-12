function [I0,Lambda] = get_Lambda_I0(B,Nph)
% Nph: Prediction horizon

nu = size(B,2);
term = diag(ones(Nph,1),0) - diag(ones(Nph-1,1),-1);
Lambda = kron(term,eye(nu));
I0 = zeros(Nph*nu,nu);
I0(1:nu,1:nu) = eye(nu);

end