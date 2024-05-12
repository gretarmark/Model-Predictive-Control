function [Phi] = get_Phi(Ad,Cdz,N)

Phi = [];
for i=1:N
Phi = [Phi ; Cdz*Ad^N];
end

end