function [Hdu,gdu,rhodu] = get_ObjectiveTermsPhidu(Nph,Qe,Lambda,I0,u_stochastic)
% Nph: Prediction horizon
% Qe:  Variance matrix

Qdu = Qe;
Wdu = chol(Qdu,'lower');
I_N = eye(Nph);
Wdu_bar = kron(I_N,Wdu);
Mdu = -(Wdu_bar*Lambda)'*Wdu_bar*I0; % Add I0 in the end
Hdu = (Wdu_bar*Lambda)'*(Wdu_bar*Lambda);  
gdu = Mdu*u_stochastic;
rhodu = (1/2)*(Wdu_bar*I0*u_stochastic)'*(Wdu_bar*I0*u_stochastic);

end