%% Research code by Agus Hasan

% Disclaimer: the analysis and results are strictly only for educational
% and research purposes and may contain errors.

% The discrete-time augmented SIRD model is used to calculate Rt. In this
% case, Rt is considered as the fifth state.

% The effective (real-time) reproduction number (Rt) is calculated based on
% extended Kalman filter. I added a low-pass filter to reduce short term
% data fluctuation.

% The contact index (CI) is calculated as Rt/Rmax

clear;
clc;

%% load data
load DATA.txt; % load data: date | month | susceptible | active cases | cummilative recovered | cummulative death

%% Infectious time
Tinf = 9;           % COVID-19 infectious time is 9 days with standard deviation of 1 day
                    % source: https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30196-1/fulltext
%%
tp  = 30;                                    % prediction time
tf  = length(DATA);                          % simulation time
N   = sum(DATA(1,3:end));                    % number of population
CFR = DATA(end,end)/(sum(DATA(end,4:6)));    % case fatality rate
td  = datetime(2020,DATA(1,2),DATA(1,1)-1) + caldays(1:tf);
tdp = datetime(2020,DATA(1,2),DATA(1,1)-1) + caldays(1:tf+tp);
dt  = 0.01;
t   = dt:dt:tf;

%% Data matrix
C = [1 0 0 0 0;     % We have data of S, I, R, and D. The fifth state,
     0 1 0 0 0;     % which is Rt, is estimated.
     0 0 1 0 0;
     0 0 0 1 0];
%% Parameters
sigma  = 1.96; %95 Confident Interval for infectious time

%% Noise
QF = diag([10 10 10 10 0.2]);   % process and measurement covariance matrices
RF = diag([100 10 10 1]);       % are considered as tuning parameters

%% For plotting
% Adding a low pass-filter to handle short-term data fluctuations 
windowSize = 500; 
b          = (1/windowSize)*ones(1,windowSize);
a          = 1;
% Plotting Rt below 1
curve11 = 0*ones(1,tf);
curve22 = 1*ones(1,tf);
x2      = [td, fliplr(td)];

%% Simulation
for j = 1:3
% Infection time
Ti     = Tinf-sigma+(j-1)*sigma;    % infection time with standard dev. 1 day
gamma  = (1-CFR)*(1/Ti);            % recovery rate
kappa  = CFR*1/Ti;                  % death rate

%% Initialization
xhat     = [N-1; 1; 0; 0; 0];   % initial condition
Pplus    = 0*eye(5);            % since we know excatly the initial conditions
xhatEff  = 0;
% for plotting
xArray       = [];
xhatArray    = [];
xhatEffArray = [];
% extended Kalman filter
for i=1:((tf-1)/dt)
     xhatArray    = [xhatArray xhat]; 
     xhatEffArray = [xhatEffArray xhatEff];      
     % assimilating the reported data
     y = [interp1(0:1:tf-1,DATA(:,3),t,'makima');
         interp1(0:1:tf-1,DATA(:,4),t,'makima');
         interp1(0:1:tf-1,DATA(:,5),t,'makima');
         interp1(0:1:tf-1,DATA(:,6),t,'makima')];
     % prediction
     xhat(1) = xhat(1)-(gamma+kappa)*xhat(5)*xhat(1)*xhat(2)*dt/N;
     xhat(2) = xhat(2)+(gamma+kappa)*xhat(5)*xhat(1)*xhat(2)*dt/N-(gamma+kappa)*xhat(2)*dt;
     xhat(3) = xhat(3)+gamma*xhat(2)*dt;
     xhat(4) = xhat(4)+kappa*xhat(2)*dt;
     xhat(5) = xhat(5);
    % calculating the Jacobian matrix
    FX    = [1-(gamma+kappa)*xhat(5)*xhat(2)*dt/N -(gamma+kappa)*xhat(5)*xhat(1)*dt/N 0 0 -(gamma+kappa)*xhat(1)*xhat(2)*dt/N;
             (gamma+kappa)*xhat(5)*xhat(2)*dt/N 1+(gamma+kappa)*xhat(5)*xhat(1)*dt/N-(gamma+kappa)*dt 0 0 (gamma+kappa)*xhat(1)*xhat(2)*dt/N;
             0 gamma*dt 1 0 0;
             0 kappa*dt 0 1 0;
             0 0 0 0 1];     
    Pmin  = FX*Pplus*FX'+QF;
    % update 
    KF    = Pmin*C'*inv(C*Pmin*C'+RF);  % Kalman gain
    xhat  = xhat + KF*(y(:,i)-C*xhat);
    Pplus = (eye(5)-KF*C)*Pmin;
    xhat(5) = max(0,xhat(5));           % the reproduction number cannot be negative
    xhatEff = (xhat(1)/N)*xhat(5);      % calculating the effective repsoduction number
end

%% Plotting

xhatArray(5,:) = filter(b,a,xhatEffArray);

xhatSArray  = [];
xhatS       = xhatArray(1,tf);
xhatIArray  = [];
xhatI       = xhatArray(2,tf);
xhatRArray  = [];
xhatR       = xhatArray(3,tf);
xhatDArray  = [];
xhatD       = xhatArray(4,tf);
xhatRtArray = [];
xhatRt      = xhatArray(5,tf);
for i=1:tf-1
    xhatSArray  = [xhatSArray xhatS];
    xhatS       = xhatArray(1,100*i);
    xhatIArray  = [xhatIArray xhatI];
    xhatI       = xhatArray(2,100*i);
    xhatRArray  = [xhatRArray xhatR];
    xhatR       = xhatArray(3,100*i);
    xhatDArray  = [xhatDArray xhatD];
    xhatD       = xhatArray(4,100*i);
    xhatRtArray = [xhatRtArray xhatRt];
    xhatRt      = xhatArray(5,100*i);
end

xhatSArray  = [xhatSArray xhatS];
xhatIArray  = [xhatIArray xhatI];
xhatRArray  = [xhatRArray xhatR];
xhatDArray  = [xhatDArray xhatD];
xhatRtArray = [xhatRtArray xhatRt];

M(j,:) = xhatRtArray;

end

curve1      = M(1,:);
xhatRtArray = M(2,:);
curve2      = M(3,:);

% RMS
RMSS = 0;
RMSI = 0;
RMSH = 0;
RMSD = 0;

for j = 1:tf
    RMSS = RMSS + sqrt(((xhatSArray(j)-DATA(j,3))/max(1,DATA(j,3)))^2);
    RMSI = RMSI + sqrt(((xhatIArray(j)-DATA(j,4))/max(1,DATA(j,4)))^2);
    RMSH = RMSH + sqrt(((xhatRArray(j)-DATA(j,5))/max(1,DATA(j,5)))^2);
    RMSD = RMSD + sqrt(((xhatDArray(j)-DATA(j,6))/max(1,DATA(j,6)))^2);
end
RMS  = (RMSS+RMSI+RMSH+RMSD)/tf;

%% Forecasting

mRt  = mean(xhatRtArray(end-5:end));
R0   = max(xhatRtArray);

cont = 0.1;

for m = 1:5

contact = cont*m;

xp   = [DATA(end,3); DATA(end,4); DATA(end,5); DATA(end,6)]; % initial condition
xpArray     = [];

cRt = contact*R0;
RtBArray = [];
for j=1:tp/dt
    RtB = mRt-((mRt-cRt)/(tp/dt))*j;
    RtBArray = [RtBArray RtB];
end
Rt = RtBArray;
for i=1:tp/dt
     xpArray = [xpArray xp]; 
     
     xp(1) = xp(1)-(gamma+kappa)*Rt(i)*xp(1)*xp(2)*dt/N;
     xp(2) = xp(2)+(gamma+kappa)*Rt(i)*xp(1)*xp(2)*dt/N-(gamma+kappa)*xp(2)*dt;
     xp(3) = xp(3)+gamma*xp(2)*dt;
     xp(4) = xp(4)+kappa*xp(2)*dt;
end

xpSArray  = [];
xpS       = xpArray(1,tf);
xpIArray  = [];
xpI       = xpArray(2,tf);
xpRArray  = [];
xpR       = xpArray(3,tf);
xpDArray  = [];
xpD       = xpArray(4,tf);
for i=1:tp
    xpSArray  = [xpSArray xpS];
    xpS       = xpArray(1,100*i);
    xpIArray  = [xpIArray xpI];
    xpI       = xpArray(2,100*i);
    xpRArray  = [xpRArray xpR];
    xpR       = xpArray(3,100*i);
    xpDArray  = [xpDArray xpD];
    xpD       = xpArray(4,100*i);
end

xIpredic(m,:) = [xhatIArray xpIArray];

end

% Plotting
figure(1)
subplot(3,1,1)
plot(td,xhatIArray,'LineWidth',6)
hold on
plot(td,DATA(:,4),'*r','LineWidth',6)
ylabel('Active Cases')
set(gca,'FontSize',24)
legend('Estimated Cases','Reported Cases','Location','northwest')
xlim([datetime(2020,DATA(1,2),DATA(1,1)), datetime(2020,DATA(end,2),DATA(end,1))])
grid on
grid minor
subplot(3,1,2)
plot(td,xhatRArray,'LineWidth',6)
hold on
plot(td,DATA(:,5),'*r','LineWidth',6)
ylabel('Recovered')
set(gca,'FontSize',24)
xlim([datetime(2020,DATA(1,2),DATA(1,1)), datetime(2020,DATA(end,2),DATA(end,1))])
grid on
grid minor
subplot(3,1,3)
plot(td,xhatDArray,'LineWidth',6)
hold on
plot(td,DATA(:,6)','*r','LineWidth',6)
ylabel('Death')
set(gca,'FontSize',24)
xlim([datetime(2020,DATA(1,2),DATA(1,1)), datetime(2020,DATA(end,2),DATA(end,1))])
grid on
grid minor

figure(2)
bar(td,[DATA(:,4) DATA(:,5) DATA(:,6)],'stacked')
set(gca,'color','none','FontSize',48)
xlim([datetime(2020,DATA(1,2),DATA(1,1)), datetime(2020,DATA(end,2),DATA(end,1))])
title('COVID-19 Cases')
legend('Active Cases','Recovered','Death','Location','northwest')
grid on
grid minor

figure(3)
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'c');
alpha(0.5)
hold on;
plot(td,xhatRtArray,'m','LineWidth',6)
hold on
inBetween = [curve11, fliplr(curve22)];
fill(x2, inBetween, 'g');
alpha(0.1)
hold on;
plot(td,ones(1,tf),'g','LineWidth',6)
set(gca,'color','none','FontSize',48)
text(tf-35,4,['Current Rt = ',num2str(xhatRtArray(end))],'FontSize',48)
ylim([0 6])
xlim([datetime(2020,DATA(1,2),DATA(1,1)), datetime(2020,DATA(end,2),DATA(end,1))])
legend('Confidence Interval 95%')
title('Real-Time Reproduction Number (Rt)')
grid on
grid minor

figure(4)
XRT = [1-(xhatRt/R0) xhatRt(end)/R0];
explode=[1 0];
labels = {'',['Current CI = ' num2str(round((xhatRt(end)/R0)*100)) '%']};
h = pie(XRT,explode,labels);
set(findobj(h,'type','text'),'fontsize',48)
title('Contact Index (CI)')
set(gca,'FontSize',48)
text(0.225,0,[datestr(datetime(2020,DATA(end,2),DATA(end,1)))],'FontSize',48,'Units','normalized','fontweight','bold')

figure(5)
plot(tdp,xIpredic(5,:),':c','LineWidth',6)
hold on;
plot(tdp,xIpredic(4,:),':g','LineWidth',6)
hold on;
plot(tdp,xIpredic(3,:),':k','LineWidth',6)
hold on;
plot(tdp,xIpredic(2,:),':b','LineWidth',6)
hold on;
plot(tdp,xIpredic(1,:),':r','LineWidth',6)
set(gca,'color','none','FontSize',48)
xline(datetime(2020,DATA(end,2),DATA(end,1)),'b','LineWidth',6)
text(tf-15,0.3*DATA(end,4),'\leftarrow Past','FontSize',48)
text(tf,0.3*DATA(end,4),'Future \rightarrow','FontSize',48)
legend_str = {'CI = 50%','CI = 40%','CI = 30%','CI = 20%','CI = 10%','Present'};
legend(legend_str(1:5),'Location','northwest')
xlim([datetime(2020,DATA(1,2),DATA(1,1)), datetime(2020,DATA(end,2),DATA(end,1)+tp)])
ylim([0 2*DATA(end,4)])
title('30-Day Forecast')
grid on
grid minor

figure(6)
sgtitle('COUNTRY','FontSize',48)
subplot(2,2,1)
bar(td,[DATA(:,4) DATA(:,5) DATA(:,6)],'stacked')
set(gca,'color','none','FontSize',24)
xlim([datetime(2020,DATA(1,2),DATA(1,1)), datetime(2020,DATA(end,2),DATA(end,1))])
title('COVID-19 Cases')
legend('Active Cases','Recovered','Death','Location','northwest')
grid on
grid minor

subplot(2,2,2)
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'c');
alpha(0.5)
hold on;
plot(td,xhatRtArray,'m','LineWidth',6)
hold on
inBetween = [curve11, fliplr(curve22)];
fill(x2, inBetween, 'g');
alpha(0.1)
hold on;
plot(td,ones(1,tf),'g','LineWidth',6)
set(gca,'color','none','FontSize',24)
text(tf-35,4,['Current Rt = ',num2str(xhatRtArray(end))],'FontSize',24)
ylim([0 6])
xlim([datetime(2020,DATA(1,2),DATA(1,1)), datetime(2020,DATA(end,2),DATA(end,1))])
legend('Confidence Interval 95%')
title('Real-Time Reproduction Number (Rt)')
grid on
grid minor

subplot(2,2,3)
XRT = [1-(xhatRt/R0) xhatRt(end)/R0];
explode=[1 0];
labels = {'',['Current CI = ' num2str(round((xhatRt(end)/R0)*100)) '%']};
h = pie(XRT,explode,labels);
set(findobj(h,'type','text'),'fontsize',24)
title('Contact Index (CI)')
set(gca,'FontSize',24)
text(0.225,-0.1,[datestr(datetime(2020,DATA(end,2),DATA(end,1)))],'FontSize',24,'Units','normalized','fontweight','bold')

subplot(2,2,4)
plot(tdp,xIpredic(5,:),':c','LineWidth',6)
hold on;
plot(tdp,xIpredic(4,:),':g','LineWidth',6)
hold on;
plot(tdp,xIpredic(3,:),':k','LineWidth',6)
hold on;
plot(tdp,xIpredic(2,:),':b','LineWidth',6)
hold on;
plot(tdp,xIpredic(1,:),':r','LineWidth',6)
set(gca,'color','none','FontSize',24)
xline(datetime(2020,DATA(end,2),DATA(end,1)),'b','LineWidth',6)
text(tf-18,0.3*DATA(end,4),'\leftarrow Past','FontSize',24)
text(tf+2,0.3*DATA(end,4),'Future \rightarrow','FontSize',24)
legend_str = {'CI = 50%','CI = 40%','CI = 30%','CI = 20%','CI = 10%','Present'};
legend(legend_str(1:5),'Location','northwest')
xlim([datetime(2020,DATA(1,2),DATA(1,1)), datetime(2020,DATA(end,2),DATA(end,1)+tp)])
ylim([0 5*DATA(end,4)])
title('30-Day Forecast')
grid on
grid minor