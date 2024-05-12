function [Qe,Rv,W,V,Qe_bar,W_bar,Gd_bar] = get_noise(N,Ed)

% Noise variance
Qe = diag([5^2 5^2]);           % Process noise
Rv = diag([2^2 2^2 2^2 2^2]);   % Measurement noise

% Generate noise vector
qe = chol(Qe,'lower');
W  = qe*randn(2,N);     % Process noise 
rv = chol(Rv,'lower');
V  = rv*randn(4,N);     % Measurement noise

% New noise W_bar
Qe_bar = diag(Ed.^2*Qe*[1; 1]);
qe_bar = chol(Qe_bar,'lower');
W_bar  = qe_bar*randn(4,N); % New process noise

Gd_bar = diag([5^2;5^2;5^2;5^2]);    % Noise matrix

end