function [t0,t_f,Ts,t,T_plot,Nsim,F1,F2,F3,F4,x0,u,d,u_rep,d_rep] = get_simulationScenario()
t0 = 0.0;           % [s] Initial time
t_f = 20*60;        % [s] Final time
Ts = 4;             % [s] Sample time
t = [t0:Ts:t_f]';   % [s] Sample instants
T_plot = t';
Nsim = length(t);

%
% t_f = 3000;
% t = [t0:Ts:t_f]';
% T_plot = t';
% Nsim = length(t);


m10 = 0.0; % [g] Liquid mass in tank 1 at time t0
m20 = 0.0; % [g] Liquid mass in tank 2 at time t0
m30 = 0.0; % [g] Liquid mass in tank 3 at time t0
m40 = 0.0; % [g] Liquid mass in tank 4 at time t0

F1 = 300; % [cm3/s] Flow rate from pump 1
F2 = 300; % [cm3/s] Flow rate from pump 2
F3 = 0.0; % Unknown stochastic variables (normal distributed) 
F4 = 0.0; % Unknown stochastic variables (normal distributed)

x0 = [m10; m20; m30; m40];

u = [F1 ; F2];
d = [F3(:,1) ; F4(:,1)];
u_rep = [repmat(F1,1,Nsim); repmat(F2,1,Nsim)];
d_rep = [repmat(F3,1,Nsim); repmat(F4,1,Nsim)];
end