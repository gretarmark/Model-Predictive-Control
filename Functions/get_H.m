function [H] = get_H(Gamma,Qz)

H = Gamma'*Qz*Gamma;

if ~issymmetric(H)
    H = (H + H')/2;
    disp('resetting H as a symmetric matrix with (H + H^T)/2')
end

end