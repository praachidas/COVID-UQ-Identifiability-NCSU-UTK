% function to solve ODEs sequentially

function f = casehosp_odes_integrateseq(P)

% initial Conditions

% SLC
y0 = zeros(1,9);ratio = [0.6 1.4];
y0(2) = P(length(P)); y0(3) = ratio(1)*P(length(P)); y0(4) = ratio(2)*P(length(P)); 
% y0(2) = 153.4551; y0(3) = ratio(1)*153.4551; y0(4) = ratio(2)*153.4551; 
y0(5) = 46; y0(6) = 130; y0(7) = 1000; y0(8) = 0; y0(9) = 0;
y0(1) = 300576 - (y0(2)+y0(3)+y0(4)+y0(5)+y0(6)+y0(7)); 
 
% % Franklin
% y0 = zeros(1,9); ratio = [0.6 1.4]; %[0.25 0.5]; 
% y0(2) = P(length(P)); y0(3) = ratio(1)*P(length(P)); y0(4) = ratio(2)*P(length(P)); 
% % y0(2) = 1.6102; y0(3) = ratio(1)*1.6102; y0(4) = ratio(2)*1.6102; %3.2999
% y0(5) = 5; y0(6) = 8; y0(7) = 64; y0(8) = 0;y0(9) = 0;
% y0(1) = 103967 - (y0(2)+y0(3)+y0(4)+y0(5)+y0(6)+y0(7));

% % SL
% y0 = zeros(1,9); ratio = [0.6 1.4];
% y0(2) = P(length(P)); y0(3) = ratio(1)*P(length(P)); y0(4) = ratio(2)*P(length(P));
% y0(5) = 214; y0(6) = 340; y0(7) = 1500; y0(8) = 0; y0(9) = 0;
% y0(1) = 994205 - (y0(2)+y0(3)+y0(4)+y0(5)+y0(6)+y0(7));

% % STC
% y0 = zeros(1,9); ratio = [0.6 1.4];
% y0(2) = P(length(P)); y0(3) = ratio(1)*P(length(P)); y0(4) = ratio(2)*P(length(P)); 
% y0(5) = 19; y0(6) = 100; y0(7) = 1000; y0(8) = 0; y0(9) = 0;
% y0(1) =  402022 - (y0(2)+y0(3)+y0(4)+y0(5)+y0(6)+y0(7));

% % Jefferson
% y0 = zeros(1,9); ratio = [0.6 1.4];
% y0(2) = P(length(P)); y0(3) = ratio(1)*P(length(P)); y0(4) = ratio(2)*P(length(P)); 
% y0(5) = 3;  y0(6) = 25; y0(7) = 500; y0(8) = 0; y0(9) = 0;
% y0(1) = 225081 - (y0(2)+y0(3)+y0(4)+y0(5)+y0(6)+y0(7));

T1 = datetime(2020,05,11):caldays(1):datetime(2020,07,09);
tvec1 = linspace(0,length(T1),length(T1))';
% opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[~,f1] = ode45(@(t,y) casehosp_odes_varybeta1(t,y,P),tvec1,y0);

% Second time period
% set initial conditions to final conditions from previous time period
y0 = zeros(1,9); y0(1) = f1(end,1); 
y0(2) = f1(end,2); y0(3) = f1(end,3); y0(4) = f1(end,4); 
y0(5) = f1(end,5); y0(6) = f1(end,6); y0(7) = f1(end,7); 
y0(8) = f1(end,8); y0(9) = f1(end,9); 

% T2 = datetime(2020,07,10):caldays(1):datetime(2020,08,19);
T2 = datetime(2020,07,09):caldays(1):datetime(2020,08,19);
tvec2 = linspace(0,length(T2),length(T2))';
% opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[~,f2] = ode45(@(t,y) casehosp_odes_varybeta2(t,y,P),tvec2,y0);

% Third time period
y0 = zeros(1,9); y0(1) = f2(end,1); 
y0(2) = f2(end,2); y0(3) = f2(end,3); y0(4) = f2(end,4); 
y0(5) = f2(end,5); y0(6) = f2(end,6); y0(7) = f2(end,7); 
y0(8) = f2(end,8); y0(9) = f2(end,9); 

% T3 = datetime(2020,08,20):caldays(1):datetime(2020,11,19);
T3 = datetime(2020,08,19):caldays(1):datetime(2020,11,19);
tvec3 = linspace(0,length(T3),length(T3))';
% opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[~,f3] = ode45(@(t,y) casehosp_odes_varybeta3(t,y,P),tvec3,y0);

% Fourth time period
y0 = zeros(1,9); y0(1) = f3(end,1); 
y0(2) = f3(end,2); y0(3) = f3(end,3); y0(4) = f3(end,4); 
y0(5) = f3(end,5); y0(6) = f3(end,6); y0(7) = f3(end,7); 
y0(8) = f3(end,8); y0(9) = f3(end,9); 

% T4 = datetime(2020,11,20):caldays(1):datetime(2020,12,27);
T4 = datetime(2020,11,19):caldays(1):datetime(2020,12,28);
tvec4 = linspace(0,length(T4),length(T4))';
% opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[~,f4] = ode45(@(t,y) casehosp_odes_varybeta4(t,y,P),tvec4,y0);

f = vertcat(f1,f2(2:end,:),f3(2:end,:),f4(2:end,:));

end