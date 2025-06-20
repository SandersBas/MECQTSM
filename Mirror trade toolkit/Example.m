%% Small number of countries
clear;clc
% Choose countries, choose from iso3_unique_sorted = {'AFG','ALB','ARE','ARG','ARM','ATG','AUS','AUT','AZE','BDI','BEL','BEN','BFA','BGD','BGR','BHS','BHR','BIH','BLR','BLZ','BOL','BRA','BRB','BRN','BWA','CAF','CAN','CHE','CHL','CHN','CIV','CMR','COG','COL','COM','CPV','CRI','CUB','CYP','CZE','DEU','DJI','DMA','DNK','DOM','DZA','ECU','EGY','ERI','ESP','EST','ETH','FIN','FJI','FRA','GAB','GBR','GEO','GHA','GIN','GMB','GNB','GRC','GRD','GTM','GUY','HND','HRV','HTI','HUN','IDN','IND','IRL','IRN','IRQ','ISL','ISR','ITA','JAM','JOR','JPN','KAZ','KEN','KGZ','KHM','KIR','KNA','KOR','KWT','LAO','LBN','LBR','LCA','LBY','LSO','LTU','LUX','LVA','MAR','MDA','MDG','MDV','MEX','MKD','MLI','MLT','MNE','MNG','MOZ','MMR','MRT','MUS','MWI','MYS','NAM','NER','NGA','NIC','NLD','NOR','NPL','NRU','NZL','OMN','PAK','PAN','PER','PHL','PLW','PNG','POL','PRK','PRT','PRY','QAT','ROU','RUS','RWA','SAU','SEN','SGP','SLB','SLE','SLV','SOM','SRB','STP','SUR','SVK','SVN','SWE','SWZ','SYC','SYR','TCD','TGO','THA','TJK','TKM','TLS','TON','TTO','TUN','TUR','TWN','TZA','UGA','UKR','URY','USA','UZB','VCT','VEN','VNM','VUT','WSM','YEM','ZAF','ZAR','ZMB','ZWE'};
countries = {'CAN','DEU','FRA','GBR','JPN','MEX','USA'};

% Choose years to use for calibration, choose from 1950-2014
years = 1994:2014;

% Choose years to produce bootstrap draws for
years_bootstrap = 2004;

[F_tilde, p, b, ME_var, P_mu, P_var] = MT_toolkit(countries,years,years_bootstrap);

%% Set of countries used in running example
clear;clc
% Choose countries, choose from iso3_unique_sorted = {'AFG','ALB','ARE','ARG','ARM','ATG','AUS','AUT','AZE','BDI','BEL','BEN','BFA','BGD','BGR','BHS','BHR','BIH','BLR','BLZ','BOL','BRA','BRB','BRN','BWA','CAF','CAN','CHE','CHL','CHN','CIV','CMR','COG','COL','COM','CPV','CRI','CUB','CYP','CZE','DEU','DJI','DMA','DNK','DOM','DZA','ECU','EGY','ERI','ESP','EST','ETH','FIN','FJI','FRA','GAB','GBR','GEO','GHA','GIN','GMB','GNB','GRC','GRD','GTM','GUY','HND','HRV','HTI','HUN','IDN','IND','IRL','IRN','IRQ','ISL','ISR','ITA','JAM','JOR','JPN','KAZ','KEN','KGZ','KHM','KIR','KNA','KOR','KWT','LAO','LBN','LBR','LCA','LBY','LSO','LTU','LUX','LVA','MAR','MDA','MDG','MDV','MEX','MKD','MLI','MLT','MNE','MNG','MOZ','MMR','MRT','MUS','MWI','MYS','NAM','NER','NGA','NIC','NLD','NOR','NPL','NRU','NZL','OMN','PAK','PAN','PER','PHL','PLW','PNG','POL','PRK','PRT','PRY','QAT','ROU','RUS','RWA','SAU','SEN','SGP','SLB','SLE','SLV','SOM','SRB','STP','SUR','SVK','SVN','SWE','SWZ','SYC','SYR','TCD','TGO','THA','TJK','TKM','TLS','TON','TTO','TUN','TUR','TWN','TZA','UGA','UKR','URY','USA','UZB','VCT','VEN','VNM','VUT','WSM','YEM','ZAF','ZAR','ZMB','ZWE'}
countries = sort({'USA', 'ARG', 'AUS', 'AUT', 'BEL', 'BEN', 'BGD', 'BOL', 'BRA', 'CAF', ...
             'CAN', 'CHE', 'CHL', 'CHN', 'CMR', 'COL', 'CRI', 'DNK', 'DOM', 'ECU', ...
             'EGY', 'ESP', 'ETH', 'FIN', 'FRA', 'GBR', 'GHA', 'GRC', 'GTM', 'HND', ...
             'IND', 'IRL', 'IRN', 'ISR', 'ITA', 'JAM', 'JOR', 'JPN', 'KEN', 'KOR', ...
             'LKA', 'MEX', 'MLI', 'MOZ', 'MUS', 'MWI', 'SGP', 'NER', 'NIC', 'NLD', ...
             'NOR', 'NPL', 'NZL', 'PAK', 'PAN', 'PER', 'PHL', 'PNG', 'PRT', 'PRY', ...
             'RWA', 'SEN', 'SLE', 'SLV', 'SWE', 'SYR', 'TGO', 'THA', 'TUN', 'TUR', ...
             'UGA', 'URY', 'VEN', 'ZAF', 'ZAR', 'ZMB', 'ZWE'});

% Choose years to use for calibration, choose from 1950-2014
years = 1986:2006;

% Choose years to produce bootstrap draws for
years_bootstrap = 1996;

[F_tilde, p, b, ME_var, P_mu, P_var] = MT_toolkit(countries,years,years_bootstrap);