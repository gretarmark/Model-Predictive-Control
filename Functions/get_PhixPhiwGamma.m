function [Phix,Phiw,Gamma] = get_PhixPhiwGamma(A,B,G,C,N)

nx = size(A,2); nu = size(B,2); nw=size(G,2); nz=size(C,1);

% Allocate memory
Nnz   = N*nz; Nnu = N*nu;
Phix  = zeros(Nnz,nx);
Phiw  = zeros(Nnz,nw);
Gamma = zeros(Nnz,Nnu);

%computer
T = C;
k1 = 1; k2=nz;
for k=1:N
    Gamma(k1:k2,1:nu) = T*B;    % C*A^(k-1)*B
    Phiw(k1:k2,1:nw)  = T*G;    % C*A^(k-1)*G   I think this is incorrect, should be C*A^(k-1)*G not C*A*G
    T = T*A;                    % C*A^k
    Phix(k1:k2,1:nx)  = T;      % C*A^k
    k1 = k1+nz; k2 = k2+nz;
end

% %Block copy
% k1 = 1; k2 = Nnz;
% k3 = 1; k4 = nu;
% for k=2:N
%     k1 = k1+nz; k2 = k2-nx;
%     k3 = k3+nu; k4 = k4+nu;
%     Gamma(k1:Nnz,k3:k4) = Gamma(1:k2,1:nu)
% end
