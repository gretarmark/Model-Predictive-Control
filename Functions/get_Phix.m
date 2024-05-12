function [Phix] = get_Phix(A,C,Nph)
% Allocate memory
nx = size(A,2);
nyz = size(C,1);
Nnz   = Nph*nyz;
Phix  = zeros(Nnz,nx);

k1 = 1; k2=nyz;
for i=1:Nph
    Phix(k1:k2,1:nx) = C*A^i;
    k1 = k1+nyz; k2 = k2+nyz;
end

end