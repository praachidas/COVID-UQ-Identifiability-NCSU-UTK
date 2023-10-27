% (Generalized) Weighted Least Squares parameter estimation and 
% uncertainty quantification for St. Louis City

clear all
close all 

tic

% Load data

data_date = (datetime(2020,05,11):caldays(7):datetime(2020,12,21))';
% Case data
data_newcases=[133;115;111;122;130;206;274;408;542;470;489;478;398;333;277;191;201;213;246;211;212;303;351;434;646;930;1154;1076;967;1040;877;882;632];
% Hospital data
data_discharges = [93;89;76;57;51;70;98;113;160;164;163;148;167;129;123;94;57;67;67;61;48;85;105;118;157;261;314;309;280;342;284;275;235];

% Create data matrix

data = [data_newcases, data_discharges];
t = data_date;
ntpts = numel(t);
T = data_date(1):caldays(1):datetime(2020,12,28);
tvec = linspace(0,length(T),length(T))';
nobvs = 2; 

% Initial estimation

k = 0; error = 1; 
%initiate first V
V = zeros(nobvs,nobvs,ntpts);
for j = 1:ntpts
    V(:,:,j) = eye(nobvs);
end

% Parameter bounds
%    beta1     beta2    beta3     beta4     kappa    rho      mu       E0
lb = [0.2      0.1      0.1       0.1       0.02     0.01     0.001    100];
ub = [0.3      0.5      0.5       0.5       0.7      0.1      0.01    250];

oldP = 0.5*(lb + ub); % set initial parameter guess
npars = length(oldP);

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

% Calculate weights
% first evaluate function (solve ODE) using estimated parameters
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
% diagonal entry being the weight for case and discharge for time point j
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

while k < 10 && error > 1e-5
    
   % update V and params for next parameter estimation
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
    [newP,f] = run(ms,problem,250); % run global search
    
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
 
    % Calculate weights
    
    % first evaluate function (solve ODE) using estimated parameters
    
    % set initial conditions
%     y0 = zeros(1,9);
%     ratio = [0.6 1.4];
%     y0(2) = newP(8);y0(3) = ratio(1)*newP(8);y0(4) = ratio(2)*newP(8); y0(5) = 46; y0(6) = 130;y0(7) = 1000;
%     y0(8) = 0; y0(9) = 0; y0(1) = 300576 - (y0(2)+y0(3)+y0(4)+y0(5)+y0(6)+y0(7)); 
%     
%     % solve ODE
%     [~,f_new] = ode45(@(t,y) casehosp_odes_varybeta(t,y,newP),tvec,y0); 
    
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

sigma1_sq = V0(1,1);
sigma2_sq = V0(2,2);
sigmas_sq3 = [sigma1_sq,sigma2_sq];

weights1_sq = zeros(ntpts,1);
weights2_sq = zeros(ntpts,1);

for j = 1:ntpts
    weights1_sq(j) = W(1,1,j);
    weights2_sq(j) = W(2,2,j);
end

P_estimated3 = newP;
solutions_estimated3 = f_new_weekly;

% Plot model and data
figure()
hold on
% plot(t,data,'o')
% plot(t,solutions_estimated3)
plot(t,data(:,1),'ro')
plot(t,data(:,2),'bo')
plot(t,solutions_estimated3(:,1),'r-')
plot(t,solutions_estimated3(:,2),'b-')
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
legend({'cases', 'discharges'},'FontSize',20,'location','northwest')
title('St. Louis City','FontSize',24)

% ------------ Residuals ------------

% (casedata(j) - model_newcases(j))/sigma1*w
% (dischargedata(j) - model_newdischarges(j))/sigma2*w
sigma1 = sqrt(sigma1_sq);
sigma2 = sqrt(sigma2_sq);
weights1_sqrt = sqrt(weights1_sq);
weights2_sqrt = sqrt(weights2_sq);

% Case data
for j = 1:ntpts
    case_weighted_residuals(j) = (data_newcases(j) - f_new_weekly(j,1))./(weights1_sqrt(j)); 
end

% Hospital data
for j = 1:ntpts
    discharges_weighted_residuals(j) = (data_discharges(j) - f_new_weekly(j,2))./(weights2_sqrt(j));  
end

% Weighted residuals vs fitted
figure()
subplot(1,2,1)
scatter(f_new_weekly(:,1), case_weighted_residuals,'LineWidth',2) %cases
hold on
plot([0 max(xlim)],[0 0],'k') %plot([0 1200],[0 0],'k')
set(gca, 'FontSize',20)
xlabel('Fitted values')
ylabel('Weighted Residuals')
title('Case Data','FontSize',22)
subplot(1,2,2)
scatter(f_new_weekly(:,2), discharges_weighted_residuals,'LineWidth',2) %discharges
set(gca, 'FontSize',20)
hold on
plot([0 max(xlim)],[0 0],'k') %plot([0 400],[0 0],'k')
xlabel('Fitted values')
ylabel('Weighted Residuals')
title('Hospital Data','FontSize',22)

% Weighted residuals vs time
figure()
subplot(1,2,1)
scatter(data_date, case_weighted_residuals,'LineWidth',2)
hold on
plot([datetime(2020,05,11) datetime(2020,12,28)],[0 0],'k')
xtickformat('MM-dd')
% xticks(data_date)
data_date_mod = (datetime(2020,05,11):caldays(21):datetime(2020,12,21))';
xticks(data_date_mod) %data_date
xtickangle(90)
set(gca,'XTickLabel',get(gca,'XTickLabel'),'FontName','Times','fontsize',18)
xlabel('Date','FontSize',22)
ylabel('Weighted Residuals','FontSize',22)
title('Case Data','FontSize',24)
subplot(1,2,2)
scatter(data_date, discharges_weighted_residuals,'LineWidth',2)
hold on
plot([datetime(2020,05,11) datetime(2020,12,28)],[0 0],'k')
xtickformat('MM-dd')
% xticks(data_date)
data_date_mod = (datetime(2020,05,11):caldays(21):datetime(2020,12,21))';
xticks(data_date_mod) %data_date
xtickangle(90)
set(gca,'XTickLabel',get(gca,'XTickLabel'),'FontName','Times','fontsize',18)
xlabel('Date','FontSize',22)
ylabel('Weighted Residuals','FontSize',22)
title('Hospital Data','FontSize',24)

%%  -------------------- Uncertainty quantification --------------------

h = 1e-10;
param = P_estimated3';
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
%     [t,x] = ode45(@casehosp_odes_varybeta,tvec,x01,[],param1);
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

% Sensitivity matrix (chi)

% State variable: confirmed cases
%chicases3_daily = dxdtheta_cases;
% for i = 1:length(param)
%     newdailycases3(:,i) = diff(chicases3_daily(:,i));
% end
chicases = dxdtheta_cases_weekly;
%FIM1 = chi1'*chi1;

% State variable: discharges
% chidischarges3_daily = dxdtheta_discharges;
% for i = 1:length(param)
%     newdailydischarges3(:,i) = diff(chidischarges3_daily(:,i));
% end
chidischarges = dxdtheta_discharges_weekly;
%FIM2 = chi2'*chi2;

% Combined sensitivity matrix
chi = zeros(ntpts,npars,nobvs);
chi(:,:,1) = chicases;
chi(:,:,2) = chidischarges;

% To find matrix M in the paper directly w/o using sensitivity matrix

% M_j = zeros(2,npars,ntpts);
% dxdtheta_M_j = zeros(2,npars,ntpts);
% 
% for k = 1:length(param) %loop through each parameter and perturb
% 
%     param1 = param;
% 
%     param1(k) = param1(k) + 1i*h;
%     [x] = casehosp_odes_integrateseq(param1);
%     weeklycases = x(find_indices,8);
%     casediff = diff(weeklycases); 
%     weeklydischarges = x(find_indices,9);
%     dischargediff = diff(weeklydischarges); 
%     
%     for l = 1:ntpts
%     dxdtheta_D_j(1,k,l) = casediff(l,:);
%     dxdtheta_D_j(2,k,l) = dischargediff(l,:);
%     end
%     
% end
% 
% M_j = imag(dxdtheta_M_j)./h;

% Plots of sensitivities
% Weekly
figure()
subplot(1,2,1)
plot(data_date,chicases);
legend({'\beta_1','\beta_2','\beta_3','\beta_4','\kappa','\rho','\mu','E_0'},'FontSize',12,'Location','northwest')
xtickformat('MM-dd')
xticks(data_date)
xtickangle(90)
set(gca,'XTickLabel',get(gca,'XTickLabel'),'FontSize',10)
xlabel('Week of','FontSize', 14)
ylabel('Parameter Sensitivity (cases)', 'FontSize', 14)
subplot(1,2,2)
plot(data_date,chidischarges);
legend({'\beta_1','\beta_2','\beta_3','\beta_4','\kappa','\rho','\mu','E_0'},'FontSize',12,'Location','northwest')
xtickformat('MM-dd')
xticks(data_date)
xtickangle(90)
set(gca,'XTickLabel',get(gca,'XTickLabel'),'FontSize',10)
xlabel('Week of','FontSize', 14)
ylabel('Parameter Sensitivity (discharges)', 'FontSize', 14)

% Covariance matrix for multiple types of data
% COV = inv(sum_j(Mj'*inv(V_hatj)*Mj))
% Each time point j will have different M and V_hat

% Define M
% First row of Mj will be partial f1 (confirmed cases) w.r.t to params 1 through p at time j
% Second row of Mj will be partial f2 (discharges) w.r.t to params 1 through p at time j
M = zeros(nobvs,npars,ntpts);
for j = 1:ntpts
    for i = 1:nobvs
        M(i,:,j) = chi(j,:,i);
    end
end

% Define V_hat
% V_hatj(11) will be V_estimate(1,1,j)
% V_hat(22) will be V_estimate(2,2,j)
V_hat = V_estimate;

% Define "inner" matrix for each j
Cov_j = zeros(npars,npars,ntpts);
for j = 1:ntpts
    Cov_j(:,:,j) = M(:,:,j)'*(inv(V_hat(:,:,j)))*M(:,:,j);
    
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
std_errors = sqrt(diag(Cov));

% Pearson correlation coefficients
[~, corr] = cov2corr(Cov);
corr = round(corr,2); % round each entry to 2 significant figures

toc

% Define cost function 
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