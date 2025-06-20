function [gamma_tilde, gamma_post_draws] = GA_toolkit(B,D_tilde)
% Input: 
%   B: number of posterior draws to produce
%   D_tilde: vector of noisy data
% Output: 
%   gamma_tilde: point estimate for gamma
%   gamma_post_draws: vector of posterior draws for gamma
% Functions needed:
%   pi_ME: function that maps D_tilde to D
%   estim: function that maps D to theta_tilde and Sigma_tilde
%   g: function that maps D and theta to gamma

%% Find point estimate for gamma
[theta_tilde, Sigma_tilde] = estim(D_tilde);
gamma_tilde = g(D_tilde,theta_tilde);

%% Generate posterior draws
gamma_post_draws = zeros(B,1);
for b = 1:B
    D_b = pi_ME(D_tilde);
    [theta_b, Sigma_tilde_b] = estim(D_b);
    theta_b_b = normrnd(theta_b, sqrt(Sigma_tilde_b),1,1);
    gamma_post_draws(b) = g(D_b, theta_b_b);
end

end

