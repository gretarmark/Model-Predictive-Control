function [Gamma] = get_Gamma(N,Ad,Bd,Czd)
% N: number of impulse response matrices to find Markov Parameters

for i=1:N
    H(:,:,i) = Czd*Ad^(i-1)*Bd;
end

Hn = size(H,1);
Hm = size(H,2);

Gamma = [];

for i=1:N
    Gamma_update = [];
    for k=N:-1:i+1
        Gamma_update = [zeros(Hn,Hm) Gamma_update];
    end
    for j=1:i
        Gamma_update = [H(:,:,j) Gamma_update];
    end
    Gamma = [Gamma ; Gamma_update];
end

end