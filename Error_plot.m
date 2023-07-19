%轨道递推
clear;
clc;

LEO_LVIM_t=load("EGM//LEO_LVIM_t.txt");LEO_LVIM_t=LEO_LVIM_t(1:400:end,1);
LEO_LVIM_X=load("EGM//LEO_LVIM_X.txt");LEO_LVIM_X=LEO_LVIM_X(1:400:end,1);
LEO_LVIM_Y=load("EGM//LEO_LVIM_Y.txt");LEO_LVIM_Y=LEO_LVIM_Y(1:400:end,1);
LEO_LVIM_Z=load("EGM//LEO_LVIM_Z.txt");LEO_LVIM_Z=LEO_LVIM_Z(1:400:end,1);

AS_LEO_LVIM_t=load("EGM//AS_LEO_LVIM_t.txt");AS_LEO_LVIM_t=AS_LEO_LVIM_t(1:100:end,1);
AS_LEO_LVIM_X=load("EGM//AS_LEO_LVIM_X.txt");AS_LEO_LVIM_X=AS_LEO_LVIM_X(1:100:end,1);
AS_LEO_LVIM_Y=load("EGM//AS_LEO_LVIM_Y.txt");AS_LEO_LVIM_Y=AS_LEO_LVIM_Y(1:100:end,1);
AS_LEO_LVIM_Z=load("EGM//AS_LEO_LVIM_Z.txt");AS_LEO_LVIM_Z=AS_LEO_LVIM_Z(1:100:end,1);

ode45_LEO=importdata("EGM//ode45_LEO.mat");
LEO_ODE_ydata=deval(ode45_LEO,LEO_LVIM_t);
LEO_sol=importdata("EGM//LEO_sol.mat");
LEO_LVIM_ydata=deval(LEO_sol,LEO_LVIM_t);
AS_LEO_LVIM_ydata=deval(LEO_sol,AS_LEO_LVIM_t);

HEO_LVIM_t=load("EGM//HEO_LVIM_t.txt");HEO_LVIM_t=HEO_LVIM_t(1:400:end,1);
HEO_LVIM_X=load("EGM//HEO_LVIM_X.txt");HEO_LVIM_X=HEO_LVIM_X(1:400:end,1);
HEO_LVIM_Y=load("EGM//HEO_LVIM_Y.txt");HEO_LVIM_Y=HEO_LVIM_Y(1:400:end,1);
HEO_LVIM_Z=load("EGM//HEO_LVIM_Z.txt");HEO_LVIM_Z=HEO_LVIM_Z(1:400:end,1);

AS_HEO_LVIM_t=load("EGM//AS_HEO_LVIM_t.txt");AS_HEO_LVIM_t=AS_HEO_LVIM_t(1:100:end,1);
AS_HEO_LVIM_X=load("EGM//AS_HEO_LVIM_X.txt");AS_HEO_LVIM_X=AS_HEO_LVIM_X(1:100:end,1);
AS_HEO_LVIM_Y=load("EGM//AS_HEO_LVIM_Y.txt");AS_HEO_LVIM_Y=AS_HEO_LVIM_Y(1:100:end,1);
AS_HEO_LVIM_Z=load("EGM//AS_HEO_LVIM_Z.txt");AS_HEO_LVIM_Z=AS_HEO_LVIM_Z(1:100:end,1);

ode45_HEO=importdata("EGM//ode45_HEO.mat");
HEO_ODE_ydata=deval(ode45_HEO,HEO_LVIM_t);
HEO_sol=importdata("EGM//HEO_sol.mat");
HEO_LVIM_ydata=deval(HEO_sol,HEO_LVIM_t);
AS_HEO_LVIM_ydata=deval(HEO_sol,AS_HEO_LVIM_t);


GEO_LVIM_t=load("EGM//GEO_LVIM_t.txt");GEO_LVIM_t=GEO_LVIM_t(1:400:end,1);
GEO_LVIM_X=load("EGM//GEO_LVIM_X.txt");GEO_LVIM_X=GEO_LVIM_X(1:400:end,1);
GEO_LVIM_Y=load("EGM//GEO_LVIM_Y.txt");GEO_LVIM_Y=GEO_LVIM_Y(1:400:end,1);
GEO_LVIM_Z=load("EGM//GEO_LVIM_Z.txt");GEO_LVIM_Z=GEO_LVIM_Z(1:400:end,1);

AS_GEO_LVIM_t=load("EGM//AS_GEO_LVIM_t.txt");AS_GEO_LVIM_t=AS_GEO_LVIM_t(1:50:end,1);
AS_GEO_LVIM_X=load("EGM//AS_GEO_LVIM_X.txt");AS_GEO_LVIM_X=AS_GEO_LVIM_X(1:50:end,1);
AS_GEO_LVIM_Y=load("EGM//AS_GEO_LVIM_Y.txt");AS_GEO_LVIM_Y=AS_GEO_LVIM_Y(1:50:end,1);
AS_GEO_LVIM_Z=load("EGM//AS_GEO_LVIM_Z.txt");AS_GEO_LVIM_Z=AS_GEO_LVIM_Z(1:50:end,1);

ode45_GEO=importdata("EGM//ode45_GEO.mat");
GEO_ODE_ydata=deval(ode45_GEO,GEO_LVIM_t);
GEO_sol=importdata("EGM//GEO_sol.mat");
GEO_LVIM_ydata=deval(GEO_sol,GEO_LVIM_t);
AS_GEO_LVIM_ydata=deval(GEO_sol,AS_GEO_LVIM_t);


figure('color',[1 1 1]);figure(1);
subplot(1, 3, 1);
semilogy(LEO_LVIM_t,sqrt((LEO_ODE_ydata(1,1:end)'-LEO_LVIM_ydata(1,1:end)').^2+(LEO_ODE_ydata(2,1:end)'-LEO_LVIM_ydata(2,1:end)').^2+(LEO_ODE_ydata(3,1:end)'-LEO_LVIM_ydata(3,1:end)').^2),'-r+','linewidth',1.5)
hold on;
semilogy(LEO_LVIM_t,sqrt((LEO_LVIM_X-LEO_LVIM_ydata(1,1:end)').^2+(LEO_LVIM_Y-LEO_LVIM_ydata(2,1:end)').^2+(LEO_LVIM_Z-LEO_LVIM_ydata(3,1:end)').^2),'-go','linewidth',1.5)
hold on;
semilogy(AS_LEO_LVIM_t,sqrt((AS_LEO_LVIM_X-AS_LEO_LVIM_ydata(1,1:end)').^2+(AS_LEO_LVIM_Y-AS_LEO_LVIM_ydata(2,1:end)').^2+(AS_LEO_LVIM_Z-AS_LEO_LVIM_ydata(3,1:end)').^2),'-b*','linewidth',1.5)
legend("ODE45","P-FAPI","AP-FAPI",'fontname','Times New Roman', 'FontSize',12,'location','SouthEast','box','off');
xlabel("Time/s",'fontname','Times New Roman');
ylabel("Position Error of LEO/m",'fontname','Times New Roman');
set(gca, 'FontSize',12); 

% annotation('arrow',[0.25 0.15],[0.8 0.75]);
% text(0.025,0.015,'\fontsize{12}\fontname{宋体}计算耗时\fontname{Times New Roman}26.434s');
% annotation('arrow',[0.25 0.15],[0.4 0.75]);
% text(0.05,0.015,'\fontsize{12}\fontname{宋体}计算耗时\fontname{Times New Roman}13.886s');

% figure('color',[1 1 1]);figure(2);
subplot(1, 3, 2);
semilogy(HEO_LVIM_t,sqrt((HEO_ODE_ydata(1,1:end)'-HEO_LVIM_ydata(1,1:end)').^2+(HEO_ODE_ydata(2,1:end)'-HEO_LVIM_ydata(2,1:end)').^2+(HEO_ODE_ydata(3,1:end)'-HEO_LVIM_ydata(3,1:end)').^2),'-r+','linewidth',1.5)
hold on;
semilogy(HEO_LVIM_t,sqrt((HEO_LVIM_X-HEO_LVIM_ydata(1,1:end)').^2+(HEO_LVIM_Y-HEO_LVIM_ydata(2,1:end)').^2+(HEO_LVIM_Z-HEO_LVIM_ydata(3,1:end)').^2),'-go','linewidth',1.5)
hold on;
semilogy(AS_HEO_LVIM_t,sqrt((AS_HEO_LVIM_X-AS_HEO_LVIM_ydata(1,1:end)').^2+(AS_HEO_LVIM_Y-AS_HEO_LVIM_ydata(2,1:end)').^2+(AS_HEO_LVIM_Z-AS_HEO_LVIM_ydata(3,1:end)').^2),'-b*','linewidth',1.5)
legend("ODE45","P-FAPI","AP-FAPI",'fontname','Times New Roman', 'FontSize',12,'location','SouthEast','box','off');
xlabel("Time/s",'fontname','Times New Roman');
ylabel("Position Error of HEO/m",'fontname','Times New Roman');
set(gca, 'FontSize',12); 


% figure('color',[1 1 1]);figure(3);
subplot(1, 3, 3);
semilogy(GEO_LVIM_t,sqrt((GEO_ODE_ydata(1,1:end)'-GEO_LVIM_ydata(1,1:end)').^2+(GEO_ODE_ydata(2,1:end)'-GEO_LVIM_ydata(2,1:end)').^2+(GEO_ODE_ydata(3,1:end)'-GEO_LVIM_ydata(3,1:end)').^2),'-r+','linewidth',1.5)
hold on;
semilogy(GEO_LVIM_t,sqrt((GEO_LVIM_X-GEO_LVIM_ydata(1,1:end)').^2+(GEO_LVIM_Y-GEO_LVIM_ydata(2,1:end)').^2+(GEO_LVIM_Z-GEO_LVIM_ydata(3,1:end)').^2),'-go','linewidth',1.5)
hold on;
semilogy(AS_GEO_LVIM_t,sqrt((AS_GEO_LVIM_X-AS_GEO_LVIM_ydata(1,1:end)').^2+(AS_GEO_LVIM_Y-AS_GEO_LVIM_ydata(2,1:end)').^2+(AS_GEO_LVIM_Z-AS_GEO_LVIM_ydata(3,1:end)').^2),'-b*','linewidth',1.5)
legend("ODE45","P-FAPI","AP-FAPI",'fontname','Times New Roman', 'FontSize',12,'location','SouthEast','box','off');
xlabel("Time/s",'fontname','Times New Roman');
ylabel("Position Error of GEO/m",'fontname','Times New Roman');
set(gca, 'FontSize',12); 
