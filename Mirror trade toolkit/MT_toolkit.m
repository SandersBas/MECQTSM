function [F_tilde, p, b, ME_var, P_mu, P_var] = MT_toolkit(countries, years, years_bootstrap)
% Input:
%   countries: choose from iso3_unique_sorted = {'AFG','ALB','ARE','ARG','ARM','ATG','AUS','AUT','AZE','BDI','BEL','BEN','BFA','BGD','BGR','BHS','BHR','BIH','BLR','BLZ','BOL','BRA','BRB','BRN','BWA','CAF','CAN','CHE','CHL','CHN','CIV','CMR','COG','COL','COM','CPV','CRI','CUB','CYP','CZE','DEU','DJI','DMA','DNK','DOM','DZA','ECU','EGY','ERI','ESP','EST','ETH','FIN','FJI','FRA','GAB','GBR','GEO','GHA','GIN','GMB','GNB','GRC','GRD','GTM','GUY','HND','HRV','HTI','HUN','IDN','IND','IRL','IRN','IRQ','ISL','ISR','ITA','JAM','JOR','JPN','KAZ','KEN','KGZ','KHM','KIR','KNA','KOR','KWT','LAO','LBN','LBR','LCA','LBY','LSO','LTU','LUX','LVA','MAR','MDA','MDG','MDV','MEX','MKD','MLI','MLT','MNE','MNG','MOZ','MMR','MRT','MUS','MWI','MYS','NAM','NER','NGA','NIC','NLD','NOR','NPL','NRU','NZL','OMN','PAK','PAN','PER','PHL','PLW','PNG','POL','PRK','PRT','PRY','QAT','ROU','RUS','RWA','SAU','SEN','SGP','SLB','SLE','SLV','SOM','SRB','STP','SUR','SVK','SVN','SWE','SWZ','SYC','SYR','TCD','TGO','THA','TJK','TKM','TLS','TON','TTO','TUN','TUR','TWN','TZA','UGA','UKR','URY','USA','UZB','VCT','VEN','VNM','VUT','WSM','YEM','ZAF','ZAR','ZMB','ZWE'};
%   years: years to use for calibration, choose from 1950-2014
%   years_bootstrap: years to produce bootstrap draws for, must be a subset of years
% Output:
%   F_tilde: vector of noisy data
%   p: vector of probabilities of true zeros
%   b: vector of probabilities of spurious zeros
%   ME_var: vector of measurement error variances
%   P_mu: vector of prior means
%   P_var: vector of prior variances
%   report: adjusted R-squared of the gravity model for the last year, and plot of log 
%           flows against log distance after partitioning out fixed effects for the last year


%% Preprocessing
% Read in mirror trade data
df = readtable('mirrortrade_dyadic_v1.csv');

% Key columns:
% -flowBAA: Bilateral trade flow from country B to country A, according to country A’s records; in million USD
% -flowBAB: Bilateral trade flow from country B to country A, according to country B’s records; in million USD

% Replace country names by iso3 codes
name_to_iso3 = containers.Map( ...
    {'Afghanistan','Albania','Algeria','Angola','Antigua & Barbuda','Argentina','Armenia','Australia','Austria','Azerbaijan',...
     'Bahamas','Bahrain','Bangladesh','Barbados','Belarus','Belgium','Belize','Benin','Bhutan','Bolivia','Bosnia and Herzegovina',...
     'Botswana','Brazil','Brunei','Bulgaria','Burkina Faso','Burundi','Cambodia','Cameroon','Canada','Cape Verde',...
     'Central African Republic','Chad','Chile','China','Colombia','Comoros','Congo','Costa Rica','Croatia','Cuba','Cyprus',...
     'Czech Republic','Czechoslovakia','Democratic Republic of the Congo','Denmark','Djibouti','Dominica','Dominican Republic',...
     'East Timor','Ecuador','Egypt','El Salvador','Equatorial Guinea','Eritrea','Estonia','Ethiopia','Fiji','Finland','France',...
     'Gabon','Gambia','Georgia','German Democratic Republic','Germany','Ghana','Greece','Grenada','Guatemala','Guinea','Guinea-Bissau',...
     'Guyana','Haiti','Honduras','Hungary','Iceland','India','Indonesia','Iran','Iraq','Ireland','Israel','Italy','Ivory Coast',...
     'Jamaica','Japan','Jordan','Kazakhstan','Kenya','Kiribati','Kuwait','Kyrgyzstan','Laos','Latvia','Lebanon','Lesotho','Liberia',...
     'Libya','Lithuania','Luxembourg','Macedonia','Madagascar','Malawi','Malaysia','Maldives','Mali','Malta','Mauritania','Mauritius',...
     'Mexico','Moldova','Mongolia','Montenegro','Morocco','Mozambique','Myanmar','Namibia','Nauru','Nepal','Netherlands','New Zealand',...
     'Nicaragua','Niger','Nigeria','North Korea','Norway','Oman','Pakistan','Palau','Panama','Papua New Guinea','Paraguay','Peru',...
     'Philippines','Poland','Portugal','Qatar','Romania','Russia','Rwanda','Samoa','Sao Tome and Principe','Saudi Arabia','Senegal',...
     'Serbia','Serbia and Montenegro','Seychelles','Sierra Leone','Singapore','Slovakia','Slovenia','Solomon Islands','Somalia',...
     'South Africa','South Korea','Spain','Sri Lanka','St. Kitts and Nevis','St. Lucia','St. Vincent and the Grenadines','Sudan',...
     'Suriname','Swaziland','Sweden','Switzerland','Syria','Taiwan','Tajikistan','Tanzania','Thailand','Togo','Tonga',...
     'Trinidad and Tobago','Tunisia','Turkey','Turkmenistan','Uganda','Ukraine','United Arab Emirates','United Kingdom',...
     'United States of America','Uruguay','Uzbekistan','Vanuatu','Venezuela','Vietnam','Yemen','Yemen Arab Republic',...
     'Yemen People''s Republic','Yugoslavia','Zambia','Zimbabwe'}, ...
    {'AFG','ALB','DZA','AGO','ATG','ARG','ARM','AUS','AUT','AZE','BHS','BHR','BGD','BRB','BLR','BEL','BLZ','BEN','BTN','BOL','BIH','BWA','BRA','BRN','BGR','BFA','BDI','KHM','CMR','CAN','CPV','CAF','TCD','CHL','CHN','COL','COM','COG','CRI','HRV','CUB','CYP','CZE','CZE','ZAR','DNK','DJI','DMA','DOM','TLS','ECU','EGY','SLV','GNQ','ERI','EST','ETH','FJI','FIN','FRA','GAB','GMB','GEO','DEU','DEU','GHA','GRC','GRD','GTM','GIN','GNB','GUY','HTI','HND','HUN','ISL','IND','IDN','IRN','IRQ','IRL','ISR','ITA','CIV','JAM','JPN','JOR','KAZ','KEN','KIR','KWT','KGZ','LAO','LVA','LBN','LSO','LBR','LBY','LTU','LUX','MKD','MDG','MWI','MYS','MDV','MLI','MLT','MRT','MUS','MEX','MDA','MNG','MNE','MAR','MOZ','MMR','NAM','NRU','NPL','NLD','NZL','NIC','NER','NGA','PRK','NOR','OMN','PAK','PLW','PAN','PNG','PRY','PER','PHL','POL','PRT','QAT','ROU','RUS','RWA','WSM','STP','SAU','SEN','SRB','SRB','SYC','SLE','SGP','SVK','SVN','SLB','SOM','ZAF','KOR','ESP','LKA','KNA','LCA','VCT','SDN','SUR','SWZ','SWE','CHE','SYR','TWN','TJK','TZA','THA','TGO','TON','TTO','TUN','TUR','TKM','UGA','UKR','ARE','GBR','USA','URY','UZB','VUT','VEN','VNM','YEM','YEM','YEM','SRB','ZMB','ZWE'});

df.countryA_iso3 = values(name_to_iso3, df.countryA_name);
df.countryB_iso3 = values(name_to_iso3, df.countryB_name);

% Only keep relevant columns
df = df(:, {'countryB_iso3', 'countryA_iso3', 'year', 'flowBAA', 'flowBAB'});
df.Properties.VariableNames = {'exp', 'imp', 'year', 'flowBAA', 'flowBAB'};

% Only keep relevant countries
df = df(ismember(df.exp, countries) & ismember(df.imp, countries), :);

% Only keep relevant years
df = df(df.year >= years(1) & df.year <= years(end), :);

% Check for missing flows and add zeros, and check for double entries
for i = 1:length(countries)
    for j = 1:length(countries)
        for k = 1:length(years)
            o = countries{i};
            d = countries{j};
            y = years(k);

            % Find matching entries
            idx = strcmp(df.exp, o) & strcmp(df.imp, d) & df.year == y;
            df_entry = df(idx, :);

            % missing flows
            if height(df_entry) < 1 && ~strcmp(o, d)
                newRow = {o, d, y, 0, 0};
                df = [df; cell2table(newRow, 'VariableNames', df.Properties.VariableNames)];
                disp({o, d, y});
            end

            % multiple entries
            if height(df_entry) > 1 && ~strcmp(o, d)
                % Remove all matching entries
                df(idx, :) = [];
                % Add a single summed entry
                flowBAA_sum = sum(df_entry.flowBAA, 'omitnan');
                flowBAB_sum = sum(df_entry.flowBAB, 'omitnan');
                newRow = {o, d, y, flowBAA_sum, flowBAB_sum};
                df = [df; cell2table(newRow, 'VariableNames', df.Properties.VariableNames)];
                disp({o, d, y});
            end
        end
    end
end

% Remove "own-country" flows
df = df(~strcmp(df.exp, df.imp), :);

% Check: sample size should now be:
expected_rows = length(countries) * (length(countries) - 1) * length(years);
disp(expected_rows);

% If for a bilateral flow one country reported only NAs and the other only non-NAs, replace the NAs with the non-NAs
for i = 1:length(countries)
    o = countries{i};
    for j = 1:length(countries)
        d = countries{j};
        if ~strcmp(o, d)
            mask_od = strcmp(df.exp, o) & strcmp(df.imp, d);
            flowBAA_vals = df.flowBAA(mask_od);
            flowBAB_vals = df.flowBAB(mask_od);

            if sum(isnan(flowBAA_vals)) == length(years) && sum(isnan(flowBAB_vals)) == 0
                df.flowBAA(mask_od) = flowBAB_vals;
                disp([o, ' ', d])
            end
            if sum(isnan(flowBAB_vals)) == length(years) && sum(isnan(flowBAA_vals)) == 0
                df.flowBAB(mask_od) = flowBAA_vals;
                disp([o, ' ', d])
            end
        end
    end
end

% Replace NAs with zeros
df.flowBAA(isnan(df.flowBAA)) = 0;
df.flowBAB(isnan(df.flowBAB)) = 0;

% Sort according to year-exporter-importer
[~, idx] = sortrows([df.year, grp2idx(categorical(df.exp)), grp2idx(categorical(df.imp))]);
df = df(idx, :);

%% Estimate Bernoulli parameters
% p: probability of a true zero
% b: probability of a spurious zero
df.p = NaN(size(df,1),1);
df.b = NaN(size(df,1),1);

for i = 1:length(countries)
    o = countries{i};
    for j = 1:length(countries)
        d = countries{j};
        if ~strcmp(o,d)
            mask_od = strcmp(df.exp, o) & strcmp(df.imp, d);
            flowBAA_vals = df.flowBAA(mask_od);
            flowBAB_vals = df.flowBAB(mask_od);

            z2 = mean((flowBAA_vals == 0) & (flowBAB_vals == 0));
            z1 = mean(((flowBAA_vals == 0) & (flowBAB_vals > 0)) | ((flowBAA_vals > 0) & (flowBAB_vals == 0)));
            z0 = mean((flowBAA_vals > 0) & (flowBAB_vals > 0));

            if z2 == 1 && z1 == 0 && z0 == 0
                df.p(mask_od) = 1;
                df.b(mask_od) = 0;
            elseif z2 == 0 && z1 == 1 && z0 == 0
                df.p(mask_od) = 0;
                df.b(mask_od) = 0.5;
            elseif z2 == 0 && z1 == 0 && z0 == 1
                df.p(mask_od) = 0;
                df.b(mask_od) = 0;
            elseif z2 > 0 && z2 < 1 && z1 > 0 && z1 < 1 && z0 == 0
                df.p(mask_od) = z2;
                df.b(mask_od) = z1;
            elseif z2 > 0 && z2 < 1 && z1 == 0 && z0 > 0 && z0 < 1
                df.p(mask_od) = z2;
                df.b(mask_od) = 0;
            elseif z2 == 0 && z1 > 0 && z1 < 1 && z0 > 0 && z0 < 1
                df.p(mask_od) = 0;
                df.b(mask_od) = z1 / (2 - z1);
            else
                df.p(mask_od) = max(1 - (z1 + 2*z0)^2 / (4*z0), 0);
                df.b(mask_od) = z1 / (z1 + 2*z0);
            end
        end
    end
end

%% Find ME variances
% Find ME variances leveraging the fact that we have a panel with two observations per year
df.MEvar = zeros(size(df,1),1);

for i = 1:length(countries)
    o = countries{i};
    for j = 1:length(countries)
        d = countries{j};
        if ~strcmp(o,d)
            mask_od = strcmp(df.exp, o) & strcmp(df.imp, d);
            flowBAA_vals = df.flowBAA(mask_od);
            flowBAB_vals = df.flowBAB(mask_od);

            % Keep only entries with positive flows in both directions
            pos_mask = (flowBAA_vals > 0) & (flowBAB_vals > 0);
            flowBAA_pos = flowBAA_vals(pos_mask);
            flowBAB_pos = flowBAB_vals(pos_mask);

            % Use model structure to estimate ME variance
            if ~isempty(flowBAA_pos)
                val = 0.5 * mean((log(flowBAA_pos) - log(flowBAB_pos)).^2);
            else
                val = 0;
            end
            df.MEvar(mask_od) = val;
        end
    end
end

%% Find prior means
% Read in distance dataset
dist_df = readtable('dist_cepii.xls');

% Keep only relevant columns
dist_df = dist_df(:, [1, 2, 12]);  % Assuming columns: iso_o, iso_d, distcap

% Only keep relevant countries
dist_df = dist_df(ismember(dist_df.iso_o, countries) & ~strcmp(dist_df.iso_o, dist_df.iso_d), :);
dist_df = dist_df(ismember(dist_df.iso_d, countries) & ~strcmp(dist_df.iso_o, dist_df.iso_d), :);

% Sort distance dataset
[~, dist_idx] = sortrows([grp2idx(categorical(dist_df.iso_o)), grp2idx(categorical(dist_df.iso_d))]);
dist_df = dist_df(dist_idx, :);

% Append distance column to df
num_years = length(unique(df.year));
df.dist = repmat(dist_df.distcap, num_years, 1);

% Rename columns
df.Properties.VariableNames = {'exp','imp','year','flowBAA','flowBAB','p','b','MEvar','dist'};

% Find prior means using model structure
df.Pmu = nan(height(df), 1);

for y = years
    % Select data for year y and positive flowBAA
    idx_y = df.year == y & df.flowBAA > 0;
    df_y = df(idx_y, :);
    
    if ~isempty(df_y)
        % Create dummy variables for exporter and importer
        [~, ~, exp_dummy] = unique(df_y.exp);
        [~, ~, imp_dummy] = unique(df_y.imp);
        exp_dummies = dummyvar(exp_dummy);
        imp_dummies = dummyvar(imp_dummy);
        
        % Run regression
        X = [log(df_y.dist), exp_dummies, imp_dummies(:,2:end)];
        y_log = log(df_y.flowBAA);
        beta = X \ y_log;
        fitted_vals = X * beta;
        
        % Store fitted values
        df.Pmu(idx_y) = fitted_vals;
    end
end

% Impute missing values with the mean across years
for i = 1:length(countries)
    for j = 1:length(countries)
        if ~strcmp(countries{i}, countries{j})
            idx = strcmp(df.exp, countries{i}) & strcmp(df.imp, countries{j});
            pmu_vals = df.Pmu(idx);
            mean_val = mean(pmu_vals(~isnan(pmu_vals)));
            df.Pmu(idx & isnan(df.Pmu)) = mean_val;
        end
    end
end

% Replace any remaining missing values with zero
df.Pmu(isnan(df.Pmu)) = 0;

%% Gather data for gravity fit report
% Compute adjusted R^2 of last year
residuals = y_log - fitted_vals;
TSS = sum((y_log - mean(y_log)).^2);
RSS = sum(residuals.^2);
R2 = 1 - RSS / TSS;
size_X_1 = size(X, 1);
size_X_2 = size(X, 2);  
R2_adj = 1 - (1 - R2) * (size_X_1 - 1) / (size_X_1 - size_X_2);

% Find residuals of log flows and log distance after partitioning out fixed effects for last year
df_y_pos = df_y(df_y.flowBAA > 0, :);
[~, ~, exp_idx] = unique(df_y_pos.exp);
[~, ~, imp_idx] = unique(df_y_pos.imp);
exp_dummies = dummyvar(exp_idx);
imp_dummies = dummyvar(imp_idx);
X_fe = [exp_dummies, imp_dummies(:,2:end)];

y_log_flow = log(df_y_pos.flowBAA);
beta_flow = X_fe \ y_log_flow;
fitted_flow = X_fe * beta_flow;
log_flow = y_log_flow - fitted_flow;

y_log_dist = log(df_y_pos.dist);
beta_dist = X_fe \ y_log_dist;
fitted_dist = X_fe * beta_dist;
log_dist = y_log_dist - fitted_dist;


%% Find prior variances
% Find prior variances using model structure
df.Pvar = nan(height(df), 1);
for i = 1:length(countries)
    o = countries{i};
    for j = 1:length(countries)
        d = countries{j};
        if ~strcmp(o,d)
            idx = strcmp(df.exp, o) & strcmp(df.imp, d);
            log_flowBAA = log(df.flowBAA(idx));
            Pmu_vals = df.Pmu(idx);
            MEvar_vals = df.MEvar(idx);
            v = var(log_flowBAA - Pmu_vals, 'omitnan') - mean(MEvar_vals, 'omitnan');
            df.Pvar(idx) = v;
        end
    end
end

% Replace missing values and negative variances with zero
df.Pvar(isnan(df.Pvar)) = 0;
df.Pvar(df.Pvar < 0) = 0;

%% Shrink MEvar and Pvar using FE model
% Create dummy variables for exporter and importer
[~, ~, exp_dummy] = unique(df.exp);
[~, ~, imp_dummy] = unique(df.imp);
exp_dummies = dummyvar(exp_dummy);
imp_dummies = dummyvar(imp_dummy);
X_dummy = [exp_dummies, imp_dummies(:,2:end)]; 

% MEvar
y = df.MEvar;
try
    b = glmfit(X_dummy, y, 'poisson', 'constant', 'off');
    df.MEvar_shrunk = exp(X_dummy * b);
catch
    df.MEvar_shrunk = zeros(height(df),1);
end

% Pvar
y = df.Pvar;
try
    b = glmfit(X_dummy, y, 'poisson', 'constant', 'off');
    df.Pvar_shrunk = exp(X_dummy * b);
catch
    df.Pvar_shrunk = zeros(height(df),1);
end

%% Clean up and collect outputs
df = df(:, {'exp','imp','year','flowBAA','flowBAB','p','b','Pmu','MEvar_shrunk','Pvar_shrunk'});
df.Properties.VariableNames = {'exp','imp','year','flowBAA','flowBAB','p','b','Pmu','MEvar','Pvar'};

% Only keep relevant years
df = df(ismember(df.year, years_bootstrap), :);

F_tilde = table2array(df(:,4));
p = table2array(df(:,6));
b = table2array(df(:,7));
ME_var = table2array(df(:,9));
P_mu = table2array(df(:,8));
P_var = table2array(df(:,10));

%% Check gravity fit
r = ksr(log_dist,log_flow,0.5);
lm = log_dist\log_flow;

scatter(log_dist,log_flow)
hold on
plot(r.x,lm.*r.x,'Color','blue','LineWidth',3)
plot(r.x,r.f,'Color','red','LineWidth',3)
hold off
set(gca,'FontSize',20)
xlabel('Log distance', 'FontSize', 20,'Interpreter','latex');
ylabel('Log trade flow', 'FontSize', 20,'Interpreter','latex');
legend('','Linear fit', 'Nonparametric fit with Gaussian kernel, bandwidth of 0.5','FontSize',20,'Interpreter','latex','Location','northwest')
title(sprintf('Log flows against log distance after partitioning out fixed effects for the year %1.0f',years(end)) )

fprintf('The adjusted R-squared of the gravity model for the year %1.0f was %4.3f\n',years(end),R2_adj) 

end

%% Auxiliary functions
function r=ksr(x,y,h,N)
% KSR   Kernel smoothing regression
%
% r=ksr(x,y) returns the Gaussian kernel regression in structure r such that
%   r.f(r.x) = y(x) + e
% The bandwidth and number of samples are also stored in r.h and r.n
% respectively.
%
% r=ksr(x,y,h) performs the regression using the specified bandwidth, h.
%
% r=ksr(x,y,h,n) calculates the regression in n points (default n=100).
%
% Without output, ksr(x,y) or ksr(x,y,h) will display the regression plot.
%
% Algorithm
% The kernel regression is a non-parametric approach to estimate the
% conditional expectation of a random variable:
%
% E(Y|X) = f(X)
%
% where f is a non-parametric function. Based on the kernel density
% estimation, this code implements the Nadaraya-Watson kernel regression
% using the Gaussian kernel as follows:
%
% f(x) = sum(kerf((x-X)/h).*Y)/sum(kerf((x-X)/h))
%
% See also gkde, ksdensity
% Example 1: smooth curve with noise
%{
x = 1:100;
y = sin(x/10)+(x/50).^2;
yn = y + 0.2*randn(1,100);
r=ksr(x,yn);
plot(x,y,'b-',x,yn,'co',r.x,r.f,'r--','linewidth',2)
legend('true','data','regression','location','northwest');
title('Gaussian kernel regression')
%}
% Example 2: with missing data
%{
x = sort(rand(1,100)*99)+1;
y = sin(x/10)+(x/50).^2;
y(round(rand(1,20)*100)) = NaN;
yn = y + 0.2*randn(1,100);
r=ksr(x,yn);
plot(x,y,'b-',x,yn,'co',r.x,r.f,'r--','linewidth',2)
legend('true','data','regression','location','northwest');
title('Gaussian kernel regression with 20% missing data')
%}
% By Yi Cao at Cranfield University on 12 March 2008.
%
% Check input and output
error(nargchk(2,4,nargin));
error(nargoutchk(0,1,nargout));
if numel(x)~=numel(y)
    error('x and y are in different sizes.');
end
x=x(:);
y=y(:);
% clean missing or invalid data points
inv=(x~=x)|(y~=y);
x(inv)=[];
y(inv)=[];
% Default parameters
if nargin<4
    N=100;
elseif ~isscalar(N)
    error('N must be a scalar.')
end
r.n=length(x);
if nargin<3
    % optimal bandwidth suggested by Bowman and Azzalini (1997) p.31
    hx=median(abs(x-median(x)))/0.6745*(4/3/r.n)^0.2;
    hy=median(abs(y-median(y)))/0.6745*(4/3/r.n)^0.2;
    h=sqrt(hy*hx);
    if h<sqrt(eps)*N
        error('There is no enough variation in the data. Regression is meaningless.')
    end
elseif ~isscalar(h)
    error('h must be a scalar.')
end
r.h=h;
% Gaussian kernel function
kerf=@(z)exp(-z.*z/2)/sqrt(2*pi);
r.x=linspace(min(x),max(x),N);
r.f=zeros(1,N);
for k=1:N
    z=kerf((r.x(k)-x)/h);
    r.f(k)=sum(z.*y)/sum(z);
end
% Plot
if ~nargout
    plot(r.x,r.f,'r',x,y,'bo')
    ylabel('f(x)')
    xlabel('x')
    title('Kernel Smoothing Regression');
end
end

