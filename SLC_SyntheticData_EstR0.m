% Monte Carlo simulations on St. Louis City synthetic dataset to calulate 
% estimates of and uncertainties in R_0 (i.e., esitmation using multiple
% replicate data sets)
% Using (Generalized) Weighted Least Squares to estimate model parameters
% Using fmincon to minimize objective function

tic
clear all
close all

rng(420)

% create time vector
data_date = (datetime(2020,05,11):caldays(7):datetime(2020,12,21))';
t = data_date;
ntpts = numel(t);
T = data_date(1):caldays(1):datetime(2020,12,28);
tvec = linspace(0,length(T),length(T))';

% set "true" parameters for generating synthetic data
params = [0.2141    0.1124   0.213   0.1701  0.0844  0.0565  0.0008 153.4551]; 

% solve ODEs using these parameters
f = casehosp_odes_integrateseq(params);

% determine weekly case and discharge data
find_indices = ismember(T,[data_date; '2020-12-28']);
weeklycases = f(find_indices,8);
casediff = diff(weeklycases);
weeklydischarges = f(find_indices,9);
dischargediff = diff(weeklydischarges);
newcases_weekly = casediff;
newdischarges_weekly = dischargediff;
t_weekly = linspace(0,length(newcases_weekly),length(newcases_weekly))';

% set gamma calculated using difference based method
gamma_cases = 0.5; 
gamma_discharges = 0.7; 
% weighting terms for errors
sigma0_cases = 1.5; % for "smaller" noise use 10^(-3);
sigma0_discharges = 0.8; % for "smaller" noise use 10^(-3);

% start loop here

m = 1;
M = 1001; %501; 
npars = length(params);
% param_estimates = zeros(M-1,npars);
param_estimates = zeros(M-1,(npars+4));
objfun_val = zeros(M-1,1);

while m < M
    
% % -----  ------ create simulated data by adding noise -----  ------

% data = x(t,r,k,x0) + x(t,r,k,x0).^gam.*(sigma0.*randn(size(t)));

% ------------- case data -------------
data_weighted_newcases_weekly = newcases_weekly + newcases_weekly.^gamma_cases.*(sigma0_cases*randn(size(t_weekly)));

% ------------- discharge data -------------
data_weighted_newdischarges_weekly = newdischarges_weekly + newdischarges_weekly.^gamma_discharges.*(sigma0_discharges*randn(size(t_weekly)));

% --------------------- Parameter estimation ---------------------
data = [data_weighted_newcases_weekly, data_weighted_newdischarges_weekly];

% --------------------- GWLS ---------------------
lb = [0.2      0.1      0.2       0.1    0.06   0.01     0.0005   150];
ub = [0.3      0.2      0.3       0.2    0.1    0.1      0.002   155];

oldP = params + 0.00002*abs(randn(size(params)));
% oldP = [params(1) + 0.002*abs(randn(1)) params(2) + 0.002*abs(randn(1)) params(3) + 0.002*abs(randn(1))...
%     params(4) + 0.002*abs(randn(1)) params(5) + 0.002*abs(randn(1)) params(6) + 0.002*abs(randn(1))...
%     params(7) + 0.0002*abs(randn(1)) params(8) + 0.02*abs(randn(1))]; 
oldP_first = oldP;
nobvs = 2;
%ntpts = length(t_weekly);
npars = length(oldP);

k = 0; error = 1; 
 %initiate first V
V = zeros(nobvs,nobvs,ntpts);
for j = 1:ntpts
    V(:,:,j) = eye(nobvs);
end

% --------------------------- fmincon ---------------------------

opts = optimoptions(@fmincon,'Algorithm','interior-point',... %set minimizing algorithm
    'StepTolerance',1e-12,... %set smaller step size (if exitflag=2)
    'FunctionTolerance',1e-12,...
    'OptimalityTolerance',1e-12,...
    'Display','off');

[newP,J] = fmincon(@(p) costfun(p,V,data,t),oldP,[],[],[],[],lb,ub,[],opts);
newP
f_new = casehosp_odes_integrateseq(newP);

data_date = (datetime(2020,05,11):caldays(7):datetime(2020,12,21))';
find_indices = ismember(T,[data_date; '2020-12-28']);
weeklycases = f_new(find_indices,8);
casediff = diff(weeklycases);
weeklydischarges = f_new(find_indices,9);
dischargediff = diff(weeklydischarges);

f_new_weekly = [casediff, dischargediff];

W = zeros(nobvs,nobvs,ntpts);

gamma1 = 0.5;
gamma2 = 0.7;

gamma = [gamma1, gamma2];

for i = 1:nobvs
    for j = 1:ntpts
        f_weekly_j(j,i) = (f_new_weekly(j,i)')^(2*gamma(i));
    end
end

for j = 1:ntpts
    W(:,:,j) = diag((f_weekly_j(j,:)')); %2x2x33 arrary
end

% redefine V as V_estimate using estimated params and W
% set initial V_estimate
V_estimate = zeros(nobvs,nobvs,ntpts);
newV = zeros(nobvs);

for j = 1:ntpts
    y_j = data(j,:)'; %2x1
    f_weekly_j = f_new_weekly(j,:)';%2x1
    vals = (y_j - f_weekly_j)*(y_j - f_weekly_j)'; % (2x1)*(1*2)
    diag_vals = diag(diag(vals)); % makes 2x2 diagonal matrix
    vals_weights = diag_vals*(inv(W(:,:,j))); % multiplying 2 2x2 diagonal matrices
    newV = newV + vals_weights; % 2x2
end

%V0 = covariance matriz

V0 = (1/(ntpts-npars))*newV;

% calculate V using V_estimate = W*V0
for j = 1:ntpts
    V_estimate(:,:,j) = (1/(ntpts-npars))*W(:,:,j)*newV;
end

params_diff = newP - oldP;
rel_error = params_diff./oldP;
error = norm(rel_error);
error_first = error;

while k < 10 && error > 1e-6
    
    % --------------------------- fmincon ---------------------------
    V = V_estimate;
    oldP = newP;
    
    opts = optimoptions(@fmincon,'Algorithm','interior-point',... %set minimizing algorithm
        'StepTolerance',1e-12,... %set smaller step size (if exitflag=2)
        'FunctionTolerance',1e-12,...
        'OptimalityTolerance',1e-12,...
        'Display','off');

   [newP,J] = fmincon(@(p) costfun(p,V,data,t),oldP,[],[],[],[],lb,ub,[],opts);
   newP
   f_new = casehosp_odes_integrateseq(newP);
   
    % find weekly new cases and discharges
    data_date = (datetime(2020,05,11):caldays(7):datetime(2020,12,21))';
    find_indices = ismember(T,[data_date; '2020-12-28']);
    weeklycases = f_new(find_indices,8);
    casediff = diff(weeklycases); 
    weeklydischarges = f_new(find_indices,9);
    dischargediff = diff(weeklydischarges);
    
   % create data matrix of model solutions for new cases and discharges
   f_new_weekly = [casediff, dischargediff];
   
   % set initial matrix for weights 
    % for each j, W is a diagonal 2x2 matrix
    % b/c we have two types of data 
    % j = 1,...,33
    W = zeros(nobvs,nobvs,ntpts);
    
    % set gamma values (obtained using difference based method)
    gamma1 = 0.5;
    gamma2 = 0.7;
    
    gamma = [gamma1, gamma2];
    
    % do the actual weights calculation
    
    % (w_j_case)^2 = (modelcasevalue_j)^(2*gamma1)
    % (w_j_hosp)^2 = (modelhospvalue_j)^(2*gamma2)
    for i = 1:nobvs
    for j = 1:ntpts
        f_weekly_j(j,i) = (f_new_weekly(j,i)')^(2*gamma(i));
    end
    end
    % f_weekly_j is jxnobvs matrix of squared weights
    % first column has squared weights for case data 
    % second column has squared weights for dishcarge data
    
    % define weight matrix
    % entries of weight matrix j should be a diagonal matrix with each
    % diagonal entry being the weight for case and discharge for time point
    % j
    for j = 1:ntpts
        W(:,:,j) = diag((f_weekly_j(j,:)')); %2x2x33 arrary
    end
    
    % redefine V as V_estimate using estimated params and W
    % set initial V_estimate
    V_estimate = zeros(nobvs,nobvs,ntpts);
    newV = zeros(nobvs);
    
    for j = 1:ntpts
        y_j = data(j,:)'; %2x1
        f_weekly_j = f_new_weekly(j,:)';%2x1
        vals = (y_j - f_weekly_j)*(y_j - f_weekly_j)'; % (2x1)*(1*2)
        diag_vals = diag(diag(vals)); % makes 2x2 diagonal matrix
        vals_weights = diag_vals*(inv(W(:,:,j))); % multiplying 2 2x2 diagonal matrices
        newV = newV + vals_weights; % 2x2
    end

    %V0 = covariance matriz
     
    V0 = (1/(ntpts-npars))*newV;
    
     % calculate V using V_estimate = W*V0
    for j = 1:ntpts
      V_estimate(:,:,j) = (1/(ntpts-npars))*W(:,:,j)*newV;
    end
    
    params_diff = newP - oldP;
    rel_error = params_diff./oldP;
    error = norm(rel_error);

    k = k+1;  
end

% Calculate R0(1) through R0(4) and store
eps = 0.3;
gamma_1 = 0.14;
gamma_2 = 0.14;
c = 1.33;
kappa = newP(5); %0.615; % 0.0844; 
mu = newP(7); % 0.0065; % 0.0008; 

for i = 1:4
    beta = newP(i);
    R0(i) = (beta*eps)/gamma_1 + (c*beta*(1-eps))/(kappa + mu + gamma_2);
end

% param_estimates(m,:) = newP;
param_estimates(m,:) = [newP,R0(1),R0(2),R0(3),R0(4)];
objfun_val(m,1) = J;
m = m + 1;
end

%% Calculate mean, std errors and correlation coefficients

mean = sum(param_estimates)./(M-1); % formula uses M, but we use M - 1 b/c m starts at 1 so M = 1001

% calculate R_0 from data, i.e., "true" R_0 values
kappa = 0.0844; % params(5)
mu = 0.0008; % params(7)
eps = 0.3;
gamma_1 = 0.14;
gamma_2 = 0.14;
c = 1.33;

for i = 1:4
    beta = params(i);
    R0_synthetic(i) = (beta*eps)/gamma_1 + (c*beta*(1-eps))/(kappa + mu + gamma_2);
end

params_all = [params,R0_synthetic];

for i = 1:(M-1)
    prod(:,:,i) = (param_estimates(i,:) - params_all)'*(param_estimates(i,:) - params_all);
end
% prod(:,:,i) is finding the product of (obsv-mean)*(obsv-mean) for all
% combinations of parameter estimates in parameter set i

sum_all = 0;

for i = 1:(M-1)
    prod_i = prod(:,:,i);
    sum_all = sum_all + prod_i;
end

cov = sum_all./(M-2);

se = sqrt(diag(cov));
CV = se./mean';

[~, corr] = cov2corr(cov);
corr = round(corr,2); % round each entry to 2 significant figures

cov_round = round(cov,6);

toc

function J = costfun(P,V,data,t)

ntpts = numel(t);
T = datetime(2020,05,11):caldays(1):datetime(2020,12,28);
tvec = linspace(0,length(T),length(T))';

f = casehosp_odes_integrateseq(P);

data_date = (datetime(2020,05,11):caldays(7):datetime(2020,12,21))';

% weekly new
find_indices = ismember(T,[data_date; '2020-12-28']);
weeklycases = f(find_indices,8);
casediff = diff(weeklycases); 
weeklydischarges = f(find_indices,9);
dischargediff = diff(weeklydischarges);
f_weekly = [casediff, dischargediff];

%sum of squared differences
%J = 0;
LS = zeros(ntpts,1);

% weekly
for j = 1:ntpts
    y_j = data(j,:)';
    f_weekly_j = f_weekly(j,:)';
    LS(j) = (y_j - f_weekly_j)'*(inv(V(:,:,j)))*(y_j - f_weekly_j);
 
end

J = sum(LS);

end
