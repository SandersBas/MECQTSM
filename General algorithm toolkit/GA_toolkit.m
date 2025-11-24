function [gamma_tilde, gamma_post_draws] = GA_toolkit(B,theta,D_tilde)
% Input: 
%   B: number of posterior draws to produce
%   theta: structural paramater
%   D_tilde: vector of noisy data
% Output: 
%   gamma_tilde: point estimate for gamma
%   gamma_post_draws: vector of posterior draws for gamma
% Functions needed:
%   pi_post: function that maps D_tilde to D
%   g: function that maps D and theta to gamma

%% Find point estimate for gamma
gamma_tilde = g(D_tilde,theta);

%% Generate posterior draws
gamma_post_draws = zeros(B,1);
for b = 1:B
    D_b = pi_ME(D_tilde);
    gamma_post_draws(b) = g(D_b, theta);
end

end

