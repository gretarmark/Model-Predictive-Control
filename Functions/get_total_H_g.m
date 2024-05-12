function [Htot,gtot] = get_total_H_g(Hz,Hdu,gz,gdu)

% Compute the total of H and g for the QP solver
% H = Hz + Hdu
% g = gz + gdu

Htot = Hz + Hdu;
gtot = gz + gdu;


end