function bk = get_bk(Phix,x_est,Phiw,w_est)
% This function should be used within a for-loop
% for step=1:Nsim-1 ...

bk = Phix*x_est + Phiw*w_est;

end