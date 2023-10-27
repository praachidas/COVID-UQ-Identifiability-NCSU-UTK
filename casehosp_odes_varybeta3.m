% function with state equations for ode system (third time period)

function dy = casehosp_odes_varybeta3(t,y,P)

%known parameters
eps = 0.3 ;
gamma_1 = 0.14;
gamma_2 = 0.14;
gamma_3 = 0.1;
gamma_4 = 0.08; %0.0836
alpha = 0.25;
c = 1.33;
mu_H = 0.01;  %0.0136

beta1 = P(1);
beta2 = P(2);
beta3 = P(3);
beta4 = P(4);
kappa = P(5);

% rho = P(5);
% mu = P(5);

rho = P(6);
mu = P(7);
% E0 = p(8);

% rho = P(5);
% mu = P(6);

% SLC
% kappa = 0.0844;
% % rho = 0.0565; 
% mu = 0.0008;
% E0 = 153.4551;
% kappa = 0.615; 
% mu = 0.0065;
% 
% kappa = 0.06; % assuming testing fraction is 0.3, mu is 0 and gamma_1 is 0.14
% % kappa = 0.2;
% mu = 0.01;

% Franklin
% kappa = 0.4368;
% % rho = 0.3200;
% %E0 = 1.6102
% % mu = 0.05;
% mu = 0.01;

% SL
% kappa = 0.1655; % example value from real data estimation
% mu = 0.0005; % example value from real data estimation

% STC
% kappa = 0.1819; % example value from real data estimation
% mu = 0.01; % upper bound real data estimation

S = y(1); E = y(2); A = y(3); I = y(4); Q = y(5); H = y(6); R = y(7);
N_1 = S+E+A+I+R;

% lambda = b(t, beta1, beta2, beta3, beta4)*(A + c*I)/N_1 ;
lambda = beta3*(A + c*I)/N_1 ;

dy = zeros(9,1);

%dS/dt
dy(1) = -lambda*S;
%dE/dt
dy(2) = lambda*S-alpha*E;
%dA/dt
dy(3) = alpha*eps*E - gamma_1*A;
%dI/dt
dy(4) = alpha*(1-eps)*E - (mu+ kappa + gamma_2)*I;
%dQ/dt
dy(5) = kappa*I - (gamma_3 + rho + mu)*Q;
%dH/dt
dy(6) = rho*Q - (mu_H + gamma_4)*H;
%dR/dt
dy(7) = gamma_1*A + gamma_2*I + gamma_3*Q + gamma_4*H;
%number of new cases
dy(8) = kappa*I;
%number leaving hospital
dy(9) = (mu_H + gamma_4)*H;

end


% %% Beta step function to vary beta over time
% 
% % start May 11 - day 1
% % End on December 27
% 
% function beta = b(t, beta1, beta2, beta3, beta4)
% 
% if t<= 60         %May 11 - July 9     
%     beta = beta1;
% elseif t>60 && t<=101  % July 10 - Aug 19
%     beta = beta2;
% elseif t>101 && t<=193 % Aug 20 - Nov 19
%     beta= beta3;
% elseif t>193    %Nov 19 - Dec 27
%     beta=beta4;
% end
%    
% end