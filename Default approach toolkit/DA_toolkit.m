function [gamma_tilde, gamma_post_draws] = DA_toolkit(B,F_tilde,p,b,ME_var,P_mu,P_var)
% Input: 
%   B: number of posterior draws to produce
%   F_tilde: vector of noisy data
%   p: vector of probabilities of true zeros
%   b: vector of probabilities of spurious zeros
%   ME_var: vector of measurement error variances
%   P_mu: vector of prior means
%   P_var: vector of prior variances
% Output: 
%   gamma_tilde: point estimate for gamma
%   gamma_post_draws: vector of posterior draws for gamma
%   report: comparison of normalized residuals histogram with standard normal pdf
% Functions needed:
%   estim: function that maps D to theta_tilde and Sigma_tilde
%   g: function that maps D and theta to gamma

%% Find point estimate for gamma
[theta_tilde, Sigma_tilde] = estim(F_tilde);
gamma_tilde = g(F_tilde,theta_tilde);

%% Find parameters of posterior
q = p./(p+b.*(1-p));

post_mean = P_var./(P_var+ME_var) .* log(F_tilde) + ME_var./(P_var+ME_var) .* P_mu;
post_var = ((1./P_var) + (1./ME_var)).^(-1);

%% Generate posterior draws
gamma_post_draws = zeros(B,1);

for b = 1:B
    prior_draw = exp(normrnd(P_mu,sqrt(P_var)));
    Q = double(rand(length(q),1)<q);
    zero_flow_draw = Q.*0 + (1-Q).*prior_draw;    

    pos_flow_draw = exp(normrnd(post_mean,sqrt(post_var)));

    F_b = pos_flow_draw;
    F_b(F_tilde==0) =  zero_flow_draw(F_tilde==0);

    [theta_b, Sigma_tilde_b] = estim(F_b);
    theta_b_b = normrnd(theta_b, sqrt(Sigma_tilde_b),1,1);
    gamma_post_draws(b) = g(F_b, theta_b_b);
end

%% Check normality
histogram((log(F_tilde)-P_mu)./sqrt(P_var+ME_var),'Normalization','pdf')
hold on
plot(sort((log(F_tilde)-P_mu)./sqrt(P_var+ME_var)),normpdf(sort((log(F_tilde)-P_mu)./sqrt(P_var+ME_var))),"LineWidth",3)
set(gca,'FontSize',20)
legend('Model residuals','Standard normal pdf', 'FontSize',20,'Interpreter','latex','Location','northwest')
xlabel('Residuals', 'FontSize', 20, 'Interpreter','latex');
ylabel('Density', 'FontSize', 20, 'Interpreter','latex');

end
