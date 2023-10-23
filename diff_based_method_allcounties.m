clear all
% close all

%% Load data

data_date = (datetime(2020,05,11):caldays(7):datetime(2020,12,21))';

% St. Louis City 
county = 'SLC';
% Case data
data_newcases=[133;115;111;122;130;206;274;408;542;470;489;478;398;333;277;191;201;213;246;211;212;303;351;434;646;930;1154;1076;967;1040;877;882;632];
% Hospital data
data_discharges = [93;89;76;57;51;70;98;113;160;164;163;148;167;129;123;94;57;67;67;61;48;85;105;118;157;261;314;309;280;342;284;275;235];

% % St. Charles County
% county = 'STC';
% % Case data
% data_newcases = [34;38;57;62;68;66;133;273;644;675;614;543;552;497;459;471;473;541;549;521;578;782;937;1037;1573;2276;2702;2400;1990;2220;2104;1938;1345];
% % Hospital data
% data_discharges = [33;31;23;36;17;31;42;62;154;201;279;235;214;198;168;186;155;168;190;185;206;240;335;335;466;880;1170;1059;809;839;813;744;581];

% % St. Louis County
% county = 'SL';
% % Case data
% data_newcases = [383;346;267;256;455;378;680;1135;1599;1444;1560;1462;1483;1305;1290;1185;1086;996;1010;1037;888;1290;1610;1809;2681;4383;5472;5052;4071;4362;3751;3417;2622];
% % Hospital data
% data_discharges = [270;268;189;169;172;181;220;379;529;644;670;651;598;539;519;482;462;369;364;374;313;442;552;670;836;1714;2282;1972;1635;1581;1582;1449;1077];
% 
% % Jefferson County
% county = 'Jefferson';
% % Case data
% data_newcases = [25;9;50;20;36;44;38;113;198;167;231;284;275;332;346;383;396;383;423;382;257;325;318;420;765;1335;1568;1492;1157;1214;1093;947;700];
% % Hospital data
% data_discharges =  [14;11;13;9;4;11;13;34;66;65;78;108;125;149;166;173;204;157;214;172;152;205;129;138;283;652;792;593;429;514;402;344;275];

% % Franklin County
% county = 'Franklin';
% % Case data
% data_newcases = [4,3,9,11,23,19,18,32,89,75,70,77,103,123,120,153,163,184,209,236,147,233,272,337,396,594,609,641,492,571,508,420,346]';
% % Hospital data
% data_discharges = [5;8;2;7;3;6;7;16;54;40;46;31;66;80;69;121;111;100;136;130;79;138;154;158;187;367;467;440;345;320;288;258;199];

y_cases = data_newcases; % vector of case data
y_hosp = data_discharges; % vector of hospital data
lengthy = length(data_date); %number of data points

%---------------------------------------------------------
%% Pick weight value (gamma)
%---------------------------------------------------------
%gamma1 = 0.3;
gamma1 = [0:0.1:1];
% gamma2 = 0.6;
gamma2 = [0:0.1:1];
% gamma2 = [1:0.1:2];
%---------------------------------------------------------
%% Pseudo Measurement Errors
%---------------------------------------------------------
% For case data
%---------------------------------------------------------
eps_cases = zeros(lengthy,1);

eps_cases(1) = (y_cases(2)-y_cases(1))/sqrt(2); %for i = 1
for i = 2:(lengthy-1)
    eps_cases(i) = (y_cases(i-1) - 2*y_cases(i) + y_cases(i+1))/(sqrt(6)); %for i=2,...,n_(k-1)
end
eps_cases(lengthy)=(y_cases(end)-y_cases(end-1))/sqrt(2); %for i = n_k

%---------------------------------------------------------
% For hospital data
%---------------------------------------------------------
eps_hosp = zeros(lengthy,1);

eps_hosp(1) = (y_hosp(2)-y_hosp(1))/sqrt(2); %for i = 1
for i = 2:(lengthy-1)
    eps_hosp(i) = (y_hosp(i-1) - 2*y_hosp(i) + y_hosp(i+1))/(sqrt(6)); %for i=2,...,n_(k-1)
end
eps_hosp(lengthy)=(y_hosp(end)-y_hosp(end-1))/sqrt(2); %for i = n_k

%---------------------------------------------------------
%% Modified Pseudo Errors
%---------------------------------------------------------
% eta_cases = zeros(lengthy,1);
% eta_hosp = zeros(lengthy,1);
% 
% for j = 1:lengthy
%     eta_cases(j) = eps_cases(j)/((abs(y_cases(j)-eps_cases(j)))^gamma1);
%     eta_hosp(j) = eps_hosp(j)/((abs(y_hosp(j)-eps_hosp(j)))^gamma2);
% end

for i = 1:length(gamma1)
    for j = 1:lengthy
        eta_cases(j,i) = eps_cases(j)/((abs(y_cases(j)-eps_cases(j)))^gamma1(i));
    end
end


for i = 1:length(gamma2)
    for j = 1:lengthy
        eta_hosp(j,i) = eps_hosp(j)/((abs(y_hosp(j)-eps_hosp(j)))^gamma2(i));
    end
end

%% Plots

% %Plot pseudo measurement errors (eps)
% %---------------------------------------------------------
% % For case data
% %---------------------------------------------------------
% figure(1)
% scatter(data_date,eps_cases)%plot eta errors
% hold on
% h = gca;set(h,'FontSize',[16]);
% plot([datetime(2020,05,11) datetime(2020,12,28)],[0 0],'k')
% ylim([-max(abs(eps_cases)),max(abs(eps_cases))])
% ylabel('\epsilon');xlabel('t')
% %title(['\gamma = ', num2str(gamma1)])
% title([num2str(county),' Case Data'])
% 
% %---------------------------------------------------------
% % For hospital data
% %---------------------------------------------------------
% figure(2)
% scatter(data_date,eps_hosp)%plot eta errors
% hold on
% h = gca;set(h,'FontSize',[16]);
% plot([datetime(2020,05,11) datetime(2020,12,28)],[0 0],'k')
% ylim([-max(abs(eps_hosp)),max(abs(eps_hosp))])
% ylabel('\epsilon');xlabel('t')
% %title(['\gamma = ', num2str(gamma2)])
% title([num2str(county),' Hospital Data'])


%Plot modified pseudo measurement errors (eps)
%---------------------------------------------------------
% For case data
%---------------------------------------------------------
% figure(3)
% %figure('visible','on')
% plot(data_date,eta_cases,'ko')%plot eta errors
% hold on
% h = gca;set(h,'FontSize',[16]);
% plot([datetime(2020,05,11) datetime(2020,12,28)],[0 0],'k')
% %xlim([min(t) max(t)])
% ylim([-max(abs(eta_cases)),max(abs(eta_cases))])
% ylabel('\eta');xlabel('t')
% title([num2str(county),' Case Data',', \gamma = ', num2str(gamma1)])

for i = 1:5
    figure(4)
    subplot(2,3,i)
    plot(data_date,eta_cases(:,i),'ko')%plot eta errors
    hold on
    h = gca;set(h,'FontSize',[16]);
    plot([datetime(2020,05,11) datetime(2020,12,28)],[0 0],'k')
    %xlim([min(t) max(t)])
    ylim([-max(abs(eta_cases(:,i))),max(abs(eta_cases(:,i)))])
    ylabel('\eta');xlabel('t')
    title(['\gamma = ', num2str(gamma1(i))])
    %title([num2str(county),' Hospital Data',', \gamma = ', num2str(gamma2(i))])
end

for i = 6:11
    figure(5)
    subplot(2,3,i-5)
    plot(data_date,eta_cases(:,i),'ko')%plot eta errors
    hold on
    h = gca;set(h,'FontSize',[16]);
    plot([datetime(2020,05,11) datetime(2020,12,28)],[0 0],'k')
    %xlim([min(t) max(t)])
    ylim([-max(abs(eta_cases(:,i))),max(abs(eta_cases(:,i)))])
    ylabel('\eta');xlabel('t')
    title(['\gamma = ', num2str(gamma1(i))])
    %title([num2str(county),' Hospital Data',', \gamma = ', num2str(gamma2(i))])
end


%---------------------------------------------------------
% For hospital data
%---------------------------------------------------------
% figure(6)
% %figure('visible','on')
% plot(data_date,eta_hosp,'ko')%plot eta errors
% hold on
% h = gca;set(h,'FontSize',[16]);
% plot([datetime(2020,05,11) datetime(2020,12,28)],[0 0],'k')
% %xlim([min(t) max(t)])
% ylim([-max(abs(eta_hosp)),max(abs(eta_hosp))])
% ylabel('\eta');xlabel('t')
% title([num2str(county),' Hospital Data',', \gamma = ', num2str(gamma2)])


for i = 1:5
    figure(7)
    subplot(2,3,i)
    plot(data_date,eta_hosp(:,i),'ko')%plot eta errors
    hold on
    h = gca;set(h,'FontSize',[16]);
    plot([datetime(2020,05,11) datetime(2020,12,28)],[0 0],'k')
    %xlim([min(t) max(t)])
    ylim([-max(abs(eta_hosp(:,i))),max(abs(eta_hosp(:,i)))])
    ylabel('\eta');xlabel('t')
    title(['\gamma = ', num2str(gamma2(i))])
    %title([num2str(county),' Hospital Data',', \gamma = ', num2str(gamma2(i))])
end

for i = 6:11
    figure(8)
    subplot(2,3,i-5)
    plot(data_date,eta_hosp(:,i),'ko')%plot eta errors
    hold on
    h = gca;set(h,'FontSize',[16]);
    plot([datetime(2020,05,11) datetime(2020,12,28)],[0 0],'k')
    %xlim([min(t) max(t)])
    ylim([-max(abs(eta_hosp(:,i))),max(abs(eta_hosp(:,i)))])
    ylabel('\eta');xlabel('t')
    title(['\gamma = ', num2str(gamma2(i))])
    %title([num2str(county),' Hospital Data',', \gamma = ', num2str(gamma2(i))])
end

%% SLC figures for paper
% 
% ---------------- Modified pseudo-errors plotted against time ----------------

% Case data

%gamma1 = 0
figure()
plot(data_date,eta_cases(:,1),'ko')%plot eta errors
hold on
plot([datetime(2020,05,11) datetime(2020,12,28)],[0 0],'k')
xtickformat('MM-dd')
% xticks(data_date)
data_date_mod = (datetime(2020,05,11):caldays(21):datetime(2020,12,21))';
xticks(data_date_mod) 
xtickangle(90)
set(gca,'XTickLabel',get(gca,'XTickLabel'),'FontName','Times','fontsize',16)
ylim([-max(abs(eta_cases(:,1))),max(abs(eta_cases(:,1)))])
ylabel('\eta','FontSize',14);xlabel('t','FontSize',18)
title(['\xi_1 = ', num2str(gamma1(1))],'FontSize',20)

%gamma1 = 0.5
figure()
plot(data_date,eta_cases(:,6),'ko')%plot eta errors
hold on
plot([datetime(2020,05,11) datetime(2020,12,28)],[0 0],'k')
xtickformat('MM-dd')
% xticks(data_date)
data_date_mod = (datetime(2020,05,11):caldays(21):datetime(2020,12,21))';
xticks(data_date_mod) 
xtickangle(90)
set(gca,'XTickLabel',get(gca,'XTickLabel'),'FontName','Times','fontsize',16)
ylim([-max(abs(eta_cases(:,6))),max(abs(eta_cases(:,6)))])
ylabel('\eta','FontSize',14);xlabel('t','FontSize',18)
title(['\xi_1 = ', num2str(gamma1(6))],'FontSize',20)

%gamma1 = 1
figure()
plot(data_date,eta_cases(:,11),'ko')%plot eta errors
hold on
plot([datetime(2020,05,11) datetime(2020,12,28)],[0 0],'k')
xtickformat('MM-dd')
% xticks(data_date)
data_date_mod = (datetime(2020,05,11):caldays(21):datetime(2020,12,21))';
xticks(data_date_mod) 
xtickangle(90)
set(gca,'XTickLabel',get(gca,'XTickLabel'),'FontName','Times','fontsize',16)
ylim([-max(abs(eta_cases(:,11))),max(abs(eta_cases(:,11)))])
ylabel('\eta','FontSize',14);xlabel('t','FontSize',18)
title(['\xi_1 = ', num2str(gamma1(11))],'FontSize',20)


% % Hospital data
% 
% % gamma2 = 0
% figure()
% plot(data_date,eta_hosp(:,1),'ko')%plot eta errors
% hold on
% plot([datetime(2020,05,11) datetime(2020,12,28)],[0 0],'k')
% xtickformat('MM-dd')
% xticks(data_date)
% xtickangle(90)
% set(gca,'XTickLabel',get(gca,'XTickLabel'),'FontName','Times','fontsize',10)
% ylim([-max(abs(eta_hosp(:,1))),max(abs(eta_hosp(:,1)))])
% ylabel('\eta','FontSize',14);xlabel('t','FontSize',14)
% title(['\xi_2 = ', num2str(gamma2(1))],'FontSize',18)
% 
% % gamma2 = 0.7
% figure()
% plot(data_date,eta_hosp(:,8),'ko')%plot eta errors
% hold on
% plot([datetime(2020,05,11) datetime(2020,12,28)],[0 0],'k')
% xtickformat('MM-dd')
% xticks(data_date)
% xtickangle(90)
% set(gca,'XTickLabel',get(gca,'XTickLabel'),'FontName','Times','fontsize',10)
% ylim([-max(abs(eta_hosp(:,8))),max(abs(eta_hosp(:,8)))])
% ylabel('\eta','FontSize',14);xlabel('t','FontSize',14)
% title(['\xi_2 = ', num2str(gamma2(8))],'FontSize',18)
% 
% % gamma2 = 1
% figure()
% plot(data_date,eta_hosp(:,11),'ko')%plot eta errors
% hold on
% plot([datetime(2020,05,11) datetime(2020,12,28)],[0 0],'k')
% xtickformat('MM-dd')
% xticks(data_date)
% xtickangle(90)
% set(gca,'XTickLabel',get(gca,'XTickLabel'),'FontName','Times','fontsize',10)
% ylim([-max(abs(eta_hosp(:,11))),max(abs(eta_hosp(:,11)))])
% ylabel('\eta','FontSize',14);xlabel('t','FontSize',14)
% title(['\xi_2 = ', num2str(gamma2(11))],'FontSize',18)

%%
% 
% % ---------------- Modified pseudo-errors plotted against data ----------------
% 
% % ---------- Specific \xi values ----------
% 
% % Case data
% 
% %gamma1 = 0
% figure()
% plot(data_newcases,eta_cases(:,1),'ko')%plot eta errors
% hold on
% % plot([datetime(2020,05,11) datetime(2020,12,28)],[0 0],'k')
% % xtickformat('MM-dd')
% % xticks(data_date)
% % xtickangle(90)
% % set(gca,'XTickLabel',get(gca,'XTickLabel'),'FontName','Times','fontsize',10)
% ylim([-max(abs(eta_cases(:,1))),max(abs(eta_cases(:,1)))])
% ylabel('\eta','FontSize',14);xlabel('Data','FontSize',14)
% title(['\xi_1 = ', num2str(gamma1(1))],'FontSize',18)
% 
% %gamma1 = 0.5
% figure()
% plot(data_newcases,eta_cases(:,6),'ko')%plot eta errors
% hold on
% % plot([datetime(2020,05,11) datetime(2020,12,28)],[0 0],'k')
% % xtickformat('MM-dd')
% % xticks(data_date)
% % xtickangle(90)
% % set(gca,'XTickLabel',get(gca,'XTickLabel'),'FontName','Times','fontsize',10)
% ylim([-max(abs(eta_cases(:,6))),max(abs(eta_cases(:,6)))])
% ylabel('\eta','FontSize',14);xlabel('Data','FontSize',14)
% title(['\xi_1 = ', num2str(gamma1(6))],'FontSize',18)
% 
% %gamma1 = 0.5
% figure()
% plot(data_newcases,eta_cases(:,11),'ko')%plot eta errors
% hold on
% % plot([datetime(2020,05,11) datetime(2020,12,28)],[0 0],'k')
% % xtickformat('MM-dd')
% % xticks(data_date)
% % xtickangle(90)
% % set(gca,'XTickLabel',get(gca,'XTickLabel'),'FontName','Times','fontsize',10)
% ylim([-max(abs(eta_cases(:,11))),max(abs(eta_cases(:,11)))])
% ylabel('\eta','FontSize',14);xlabel('data','FontSize',14)
% title(['\xi_1 = ', num2str(gamma1(11))],'FontSize',18)

 
% ---------- Multiple \xi values ----------

% -------- Case data --------

% for i = 1:5
%     figure(9)
%     subplot(2,3,i)
%     plot(y_cases,eta_cases(:,i),'ko') %plot eta errors
%     hold on
%     plot([0 max(xlim)],[0 0],'k')
%     h = gca;set(h,'FontSize',[16]);
%     ylim([-max(abs(eta_cases(:,i))),max(abs(eta_cases(:,i)))])
%     ylabel('\eta');xlabel('New cases')
%     title(['\gamma = ', num2str(gamma1(i))])
%     %title([num2str(county),' Hospital Data',', \gamma = ', num2str(gamma2(i))])
% end
% 
% for i = 6:11
%     figure(10)
%     subplot(2,3,i-5)
%     plot(y_cases,eta_cases(:,i),'ko')%plot eta errors
%     hold on
%     plot([0 max(xlim)],[0 0],'k')
%     h = gca;set(h,'FontSize',[16]);
%     ylim([-max(abs(eta_cases(:,i))),max(abs(eta_cases(:,i)))])
%     ylabel('\eta');xlabel('New cases')
%     title(['\gamma = ', num2str(gamma1(i))])
%     %title([num2str(county),' Hospital Data',', \gamma = ', num2str(gamma2(i))])
% end

% -------- Hospital data --------

% for i = 1:5
%     figure(13)
%     subplot(2,3,i)
%     plot(y_hosp,eta_hosp(:,i),'ko')%plot eta errors
%     hold on
%     plot([0 max(xlim)],[0 0],'k')
%     h = gca;set(h,'FontSize',[16]);
%     ylim([-max(abs(eta_hosp(:,i))),max(abs(eta_hosp(:,i)))])
%     ylabel('\eta');xlabel('New discharges')
%     title(['\gamma = ', num2str(gamma2(i))])
%     %title([num2str(county),' Hospital Data',', \gamma = ', num2str(gamma2(i))])
% end
% 
% for i = 6:11
%     figure(14)
%     subplot(2,3,i-5)
%     plot(y_hosp,eta_hosp(:,i),'ko')%plot eta errors
%     hold on
%     plot([0 max(xlim)],[0 0],'k')
%     h = gca;set(h,'FontSize',[16]);
%     ylim([-max(abs(eta_hosp(:,i))),max(abs(eta_hosp(:,i)))])
%     ylabel('\eta');xlabel('New discharges')
%     title(['\gamma = ', num2str(gamma2(i))])
%     %title([num2str(county),' Hospital Data',', \gamma = ', num2str(gamma2(i))])
% end