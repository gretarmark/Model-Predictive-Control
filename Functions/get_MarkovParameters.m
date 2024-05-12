function [H] = get_MarkovParameters(A,B,C,k)
% for k=0,         H_k = 0
% for k=1,2,...    H_k = C*A^(k-1)*B

if k==0
    H(1) = 0;
else
    for i=1:k+1
        H(:,:,i) = C*A^(i-1)*B;
    end
end