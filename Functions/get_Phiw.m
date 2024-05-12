function [Phiw] = get_Phiw(A,C,G,N)

nw = size(G,2);
nyz = size(C,1);
Nnz = N*nyz;
Phiw  = zeros(Nnz,nw);

k1 = 1; k2=nyz;
for i=1:N
    Phiw(k1:k2,1:nw)  = C*A^(i-1)*G;
    k1 = k1+nyz; k2 = k2+nyz;
end

end