% (Generalized) Weighted Least Squares parameter estimation and 
% uncertainty quantification for St. Louis City using synthetic data
% Using MultiStart to minimize objective function


tic 

clear all 
close all

rng(420)

% Solve ODE with fixed parameters
data_date = (datetime(2020,05,11):caldays(7):datetime(2020,12,21))';
t = data_date;
ntpts = numel(t);
T = data_date(1):caldays(1):datetime(2020,12,28);
tvec = linspace(0,length(T),length(T))';

% Set true parameter values
%         beta1     beta2    beta3   beta4    kappa   rho     mu      E0
params = [0.2141    0.1124   0.213   0.1701   0.0844  0.0565 0.0008  153.4551];

% ----------- Integrating sequentially -----------

f = casehosp_odes_integrateseq(params);

find_indices = ismember(T,[data_date; '2020-12-28']);
weeklycases = f(find_indices,8);
casediff = diff(weeklycases);
weeklydischarges = f(find_indices,9);
dischargediff = diff(weeklydischarges);
newcases_weekly = casediff;
newdischarges_weekly = dischargediff;
t_weekly = linspace(0,length(newcases_weekly),length(newcases_weekly))';

gamma_cases = 0.5; %10^(-5)
gamma_discharges = 0.7; %10^(-5)
sigma0_cases = 1.5; % use 10^(-3) for "small" noise; 
sigma0_discharges = 0.8; % 10^(-3)

% create synthetic data set with noise 

data_weighted_newcases_weekly = newcases_weekly + (newcases_weekly.^gamma_cases).*(sigma0_cases*randn(size(t_weekly)));
data_weighted_newdischarges_weekly = newdischarges_weekly + (newdischarges_weekly.^gamma_discharges).*(sigma0_discharges*randn(size(t_weekly)));
data = [data_weighted_newcases_weekly, data_weighted_newdischarges_weekly];
nobvs = 2; 

figure()
scatter(data_date,data_weighted_newcases_weekly)
hold on
scatter(data_date, data_weighted_newdischarges_weekly)

k = 0; error = 1; 
%initiate first V
V = zeros(nobvs,nobvs,ntpts);
for j = 1:ntpts
     V(:,:,j) = eye(nobvs);
%      V(:,:,j) = V;
end

oldP = params + 0.0002*abs(randn(size(params)));
% oldP = [params(1) + 0.002*abs(randn(1)) params(2) + 0.002*abs(randn(1)) params(3) + 0.002*abs(randn(1))...
%     params(4) + 0.002*abs(randn(1)) params(5) + 0.002*abs(randn(1)) params(6) + 0.002*abs(randn(1))...
%     params(7) + 0.00002*abs(randn(1)) params(8) + 0.002*abs(randn(1))];
oldP_first = oldP; 
nobvs = 2;
ntpts = length(t_weekly);
npars = length(oldP);

lb = [0.2      0.1      0.2       0.1    0.06   0.01     0.0005   150];
ub = [0.3      0.2      0.3       0.2    0.1    0.1      0.002   155];

% --------------------------- MultiStart ---------------------------

opts = optimoptions(@fmincon,'Algorithm','interior-point',... %set minimizing algorithm
        'StepTolerance',1e-12,... %set smaller step size (if exitflag=2) 1e-6
        'FunctionTolerance',1e-12,...
        'OptimalityTolerance',1e-12); % 1e-8
%         'StepTolerance',1e-5,... %set smaller step size (if exitflag=2)
%         'ConstraintTolerance',1e-6); %set smaller constraint tolerance (if exitflag=2)
problem = createOptimProblem('fmincon',... %minimizer
    'objective',@(p) costfun3(p,V,data,t),... %optimization function
    'x0',oldP,...
    'lb',lb,...
    'ub',ub,...
    'options',opts);
ms = MultiStart(MultiStart,'Display','iter');
[newP,f] = run(ms,problem,1);

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

% % update V and params for next parameter estimation
% V = V_estimate;
% oldP = newP;

while k < 10 && error > 1e-6
    
     % -------------------- Using Multistart -------------------- 
   
    V = V_estimate;
    oldP = newP;
    
     if k < 1
    opts = optimoptions(@fmincon,'Algorithm','interior-point',... %set minimizing algorithm
        'StepTolerance',1e-12,... %set smaller step size (if exitflag=2) 1e-6
        'FunctionTolerance',1e-12,...
        'OptimalityTolerance',1e-12); % 1e-8
%         'StepTolerance',1e-5,... %set smaller step size (if exitflag=2)
%         'ConstraintTolerance',1e-6); %set smaller constraint tolerance (if exitflag=2)
    problem = createOptimProblem('fmincon',... %minimizer
        'objective',@(p) costfun3(p,V,data,t),... %optimization function
        'x0',oldP,...
        'lb',lb,...
        'ub',ub,...
        'options',opts);
    ms = MultiStart(MultiStart,'Display','iter'); 
    [newP,f] = run(ms,problem,250); % 500
    
   else
       
    opts = optimoptions(@fmincon,'Algorithm','interior-point',... %set minimizing algorithm
     'StepTolerance',1e-12,... %set smaller step size (if exitflag=2) 1e-6
        'FunctionTolerance',1e-12,...
        'OptimalityTolerance',1e-12); % 1e-8
%         'StepTolerance',1e-5,... %set smaller step size (if exitflag=2)
%         'ConstraintTolerance',1e-6); %set smaller constraint tolerance (if exitflag=2)
    problem = createOptimProblem('fmincon',... %minimizer
        'objective',@(p) costfun3(p,V,data,t),... %optimization function
        'x0',oldP,...
        'lb',lb,...
        'ub',ub,...
        'options',opts);
    ms = MultiStart(MultiStart,'Display','iter');
    [newP,f] = run(ms,problem,10); 
     end
   
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
  
    % update V and params for next parameter estimation
%     V = V_estimate;
%     oldP = newP;
    k = k+1;
end

K = k;
sigma1_sq = V0(1,1);
sigma2_sq = V0(2,2);

weights1_sq = zeros(ntpts,1);
weights2_sq = zeros(ntpts,1);

for j = 1:ntpts
    weights1_sq(j) = W(1,1,j);
    weights2_sq(j) = W(2,2,j);
end

P_estimated = newP;
solutions_estimated = f_new_weekly;

f_true = casehosp_odes_integrateseq(params);
weeklycases_true = f_true(find_indices,8);
casediff_true  = diff(weeklycases_true);
weeklydischarges_true  = f_true(find_indices,9);
dischargediff_true  = diff(weeklydischarges_true);

f_true_weekly = [casediff_true, dischargediff_true];
solutions_estimated_true = f_true_weekly;

% P_estimated = newP_fminsearch;
% modelsol_fminsearch = casehosp_odes_integrateseq(newP_fminsearch); %f_new_weekly;
%  
% find weekly new cases and discharges
% weeklycases = f_new(find_indices,8);
% casediff = diff(weeklycases);
% weeklydischarges = f_new(find_indices,9);
% dischargediff = diff(weeklydischarges);
%     
% % create data matrix of model solutions for new cases and discharges
% f_new_weekly = [casediff, dischargediff];
% solutions_estimated = f_new_weekly;

%Plot model and data
figure()
hold on
plot(t,data,'o')
% plot(t,solutions_estimated)
% plot(t,solutions_estimated_true)
plot(t,solutions_estimated(:,1),'k-')
plot(t,solutions_estimated_true(:,1),'k--')
plot(t,solutions_estimated(:,2),'k-')
plot(t,solutions_estimated_true(:,2),'k--')
plot([datetime(2020,07,09) datetime(2020,07,09)],[0 1200], 'k')
plot([datetime(2020,08,19) datetime(2020,08,19)],[0 1200], 'k')
plot([datetime(2020,11,19) datetime(2020,11,19)],[0 1200], 'k')
xtickformat('MM-dd')
%xticks(data_date)
data_date_mod = (datetime(2020,05,11):caldays(21):datetime(2020,12,21))';
xticks(data_date_mod) %data_date
xtickangle(90)
set(gca,'XTickLabel',get(gca,'XTickLabel'),'FontName','Times','fontsize',18)
xlabel('Date','FontSize',22)
ylabel('Number of individuals','FontSize',22)
% legend({'synthetic case data', 'synthetic discharge data'},'FontSize',14,'location','northwest')
legend({'synthetic case data', 'synthetic discharge data', 'estimated','true'},'FontSize',18,'location','northwest')
title('St. Louis City','FontSize',24)

% % -------------------- Uncertainty quantification --------------------
% 
h = 1e-10;
param = params; % P_estimated';
% x0 = y0;

dxdtheta_cases = zeros(length(T),npars); %initialize for storage
dxdtheta_discharges = zeros(length(T),npars); 
dxdtheta_cases_weekly = zeros(length(data_date),npars); %initialize for storage
dxdtheta_discharges_weekly = zeros(length(data_date),npars); 

for k = 1:length(param) %loop through each parameter and perturb

    param1 = param;

    param1(k) = param1(k) + 1i*h;
    
%     x01 = x0;
%     x01(2) = param1(npars); % setting E(0) = perturbed/non-perturbed value
    % [t,x] = ode45(@casehosp_odes_varybeta,tvec,x01,[],param1); 
    [x] = casehosp_odes_integrateseq(param1);
    dxdtheta_cases(:,k) = x(:,8);
    dxdtheta_discharges(:,k) = x(:,9);
    weeklycases = x(find_indices,8);
    casediff = diff(weeklycases); 
    dxdtheta_cases_weekly(:,k) = casediff;
    weeklydischarges = x(find_indices,9);
    dischargediff = diff(weeklydischarges); 
    dxdtheta_discharges_weekly(:,k) = dischargediff;
end   

%Calculate partials
dxdtheta_cases =  imag(dxdtheta_cases)./h;
dxdtheta_cases_weekly = imag(dxdtheta_cases_weekly)./h; % This is chi1
dxdtheta_discharges =  imag(dxdtheta_discharges)./h;
dxdtheta_discharges_weekly = imag(dxdtheta_discharges_weekly)./h; % This is chi2

% Sensitivity matrix

% State variable: confirmed cases
chicases = dxdtheta_cases_weekly;
FIM1 = chicases'*chicases;
cond(FIM1)

% State variable: discharges
chidischarges = dxdtheta_discharges_weekly;
FIM2 = chidischarges'*chidischarges;
cond(FIM2)

% Combined sensitivity matrix
chi = zeros(ntpts,npars,nobvs);
chi(:,:,1) = chicases;
chi(:,:,2) = chidischarges;

% Plots of sensitivities
% Weekly
figure()
subplot(1,2,1)
plot(data_date,chicases, 'LineWidth', 2);
legend({'\beta_1','\beta_2','\beta_3','\beta_4','\kappa','\rho','\mu','E_0'},'FontSize',12,'Location','northwest')
xtickformat('MM-dd')
xticks(data_date_mod) % xticks(data_date)
xtickangle(90)
set(gca,'XTickLabel',get(gca,'XTickLabel'),'FontSize',18)
xlabel('Week of','FontSize', 22)
ylabel('Parameter Sensitivity (cases)', 'FontSize', 22)
subplot(1,2,2)
plot(data_date,chidischarges, 'LineWidth', 2);
legend({'\beta_1','\beta_2','\beta_3','\beta_4','\kappa','\rho','\mu','E_0'},'FontSize',12,'Location','northwest')
xtickformat('MM-dd')
xticks(data_date_mod) % xticks(data_date)
xtickangle(90)
set(gca,'XTickLabel',get(gca,'XTickLabel'),'FontSize',18)
xlabel('Week of','FontSize', 22)
ylabel('Parameter Sensitivity (discharges)', 'FontSize', 22)

% Define D (denoted as M in the manuscript)
% First row of Dj will be partial f1 (confirmed cases) w.r.t to params 1 through p at time j
% Second row of Dj will be partial f2 (discharges) w.r.t to params 1 through p at time j
D = zeros(nobvs,npars,ntpts);
for j = 1:ntpts
    for i = 1:nobvs
        D(i,:,j) = chi(j,:,i);
    end
end

% Define V_hat
% V_hatj(11) will be V_estimate(1,1,j)
% V_hat(22) will be V_estimate(2,2,j)
V_hat = V_estimate;
% V_hat = V;

% Define "inner" matrix for each j
Cov_j = zeros(npars,npars,ntpts);
for j = 1:ntpts
    Cov_j(:,:,j) = D(:,:,j)'*(inv(V_hat(:,:,j)))*D(:,:,j);
end

% Sum values of Cov_j for j = 1:33
Cov_sum = 0;
for j = 1:ntpts
    Cov_sum_j = Cov_j(:,:,j);
    Cov_sum = Cov_sum + Cov_sum_j;
end

% Calculate covariance matrix
Cov = inv(Cov_sum);

% Calculate standard errors
std_errors = sqrt(diag(Cov))

% Calculate coefficients of variation
CV = std_errors./params' % newP';

% Pearson correlation coefficients
[~, corrfull] = cov2corr(Cov)
corr = round(corrfull,2) % round each entry to 2 significant figures

%% Structural Identifiability 
% Define D2 (look at the rank of this for structural identifiability)
% basically concatenating rows of D

splitD = num2cell(D, [1 2]); %split D keeping dimension 1 and 2 intact
D2 = vertcat(splitD{:}); %concatenate rows of D

rank(D2)
[U,S,V]=svd(D2);
SingularValues=diag(S)';
Log10SingularValues=log10(SingularValues);
RowTable.SingularValues={Log10SingularValues};
RowTable.MinSingVal=min(min(SingularValues));
minS=min(SingularValues(:,end));
maxS=max(SingularValues(:,1));

singularvalues_final=log10(SingularValues)';

figure()
plot(singularvalues_final,'.','MarkerSize',32);
xlabel('Singular values in descending order');
ylabel('Log_{10}(Singular Values)');
title('^{10}Log(Singular Values)','FontWeight','bold');

toc

function J3 = costfun3(P,V,data,t)

ntpts = numel(t);
T = datetime(2020,05,11):caldays(1):datetime(2020,12,28);
tvec = linspace(0,length(T),length(T))';

f = casehosp_odes_integrateseq(P);

data_date = (datetime(2020,05,11):caldays(7):datetime(2020,12,21))';
find_indices = ismember(T,[data_date; '2020-12-28']);
weeklycases = f(find_indices,8);
casediff = diff(weeklycases); 
weeklydischarges = f(find_indices,9);
dischargediff = diff(weeklydischarges);
f_weekly = [casediff, dischargediff];

%sum of squared differences
LS = zeros(ntpts,1);

for j = 1:ntpts
    y_j = data(j,:)';
    f_weekly_j = f_weekly(j,:)';
    LS(j) = (y_j - f_weekly_j)'*(inv(V(:,:,j)))*(y_j - f_weekly_j);
 
end

J3 = sum(LS);

end
