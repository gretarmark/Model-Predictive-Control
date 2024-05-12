function [y_stochastic,y_est] = KalmanFilterDynamic(x0,xs,ys,us,N,u_rep,Ad,Bd,Ed,C,W,d_rep,V,Qe_bar,Gd_bar,Rv)



% static KF (if you want to test you can replace dynamic kf gain with static kalman filter gain )
% P_cov = dare(Ad', C', Qe_bar, Rv)
% R_e = C*P_cov*C' + Rv;
% fprintf('Kalman gain is:')
% K_f = P_cov*(C'/R_e)

% dynamic kf
P_d = diag([0.015^2, 0.015^2, 0.015^2, 0.015^2]) % initial guess on P 

% pre-allocate memory for storing results 
x_stochastic = zeros(4,N);
x_stochastic(:,1) = x0 - xs;
xbar_stochastic = x_stochastic; % simulation with W_bar

uss = [300 ; 300];
u_stochastic = u_rep - uss;
y_stochastic = zeros(4,N) - ys; % y measurement 
y_est = zeros(4, N)- ys;

ybar_stochastic = y_stochastic; % simulation with W_bar

% simulation for normal case 
for i=1:N-1
   P_d = Ad*P_d*Ad' + Gd_bar*Qe_bar*Gd_bar';
   K_df = P_d*C'*(C*P_d*C'+Rv); % kalman gain 
   % simulation 
   x_stochastic(:,i+1) = Ad*x_stochastic(:,i) + Bd*u_stochastic(:,i) + Ed*(W(:,i)+d_rep(:,i));
   y_stochastic(:,i+1) = C*x_stochastic(:,i+1) + V(:,i+1);
   e = y_stochastic(:,i+1) - C*x_stochastic(:,i+1);
   x_est = x_stochastic(:, i+1) + K_df*e;
   y_est(:, i+1) = C*x_est;
   P_d = P_d - K_df*(C*P_d*C'+Rv) *K_df';
end



end