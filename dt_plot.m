clear;clc;
dt_LEO=load("EGM//dt_LEO.txt");
ic1_LEO=load("EGM//ic1_LEO.txt");
ic2_LEO=load("EGM//ic2_LEO.txt");
t_LEO=1:1:length(dt_LEO);


figure('color',[1 1 1]);figure(1);
subplot(1, 3, 1);
yyaxis left;
scatter(t_LEO,dt_LEO','filled','Marker','o');
ylim([200,3500]);
ylabel('$dt/s$', 'Interpreter', 'latex');
hold on;
topline=yline(3500,'-k','LineWidth', 3);
set(topline, 'HandleVisibility', 'off'); % 设置直线不参与图例的排序
yyaxis right;
scatter(t_LEO,ic1_LEO,'filled','Marker','^');
hold on;
scatter(t_LEO,ic2_LEO,'filled','Marker','v');
ylim([0,42]);
ylabel('Number of iterations','FontName','Times New Roman');
xlabel('Ordinal number of step','FontName','Times New Roman');
legend("$dt$","$I_{c1} $","$I_{c2} $", 'Interpreter', 'latex');


dt_HEO=load("EGM//dt_HEO.txt");
ic1_HEO=load("EGM//ic1_HEO.txt");
ic2_HEO=load("EGM//ic2_HEO.txt");
t_HEO=1:1:length(dt_HEO);


% figure('color',[1 1 1]);figure(2);
subplot(1, 3, 2);
yyaxis left;
scatter(t_HEO,dt_HEO','filled','Marker','o');
ylim([500,10000]);
ylabel('$dt/s$', 'Interpreter', 'latex');
hold on;
topline=yline(10000,'-k','LineWidth', 3);
set(topline, 'HandleVisibility', 'off'); % 设置直线不参与图例的排序
yyaxis right;
scatter(t_HEO,ic1_HEO,'filled','Marker','^');
hold on;
scatter(t_HEO,ic2_HEO,'filled','Marker','v');
ylim([0,42]);
ylabel('Number of iterations','FontName','Times New Roman');
xlabel('Ordinal number of step','FontName','Times New Roman');
legend("$dt$","$I_{c1} $","$I_{c2} $", 'Interpreter', 'latex');

dt_GEO=load("EGM//dt_GEO.txt");
ic1_GEO=load("EGM//ic1_GEO.txt");
ic2_GEO=load("EGM//ic2_GEO.txt");
t_GEO=1:1:length(dt_GEO);


% figure('color',[1 1 1]);figure(3);
subplot(1, 3, 3);
xlim([0,20]);
yyaxis left;
scatter(t_GEO,dt_GEO','filled','Marker','o');
%ylim([500,10000]);
ylabel('$dt/s$', 'Interpreter', 'latex');
hold on;
topline=yline(80000,'-k','LineWidth', 3);
set(topline, 'HandleVisibility', 'off'); % 设置直线不参与图例的排序
yyaxis right;
scatter(t_GEO,ic1_GEO,'filled','Marker','^');
hold on;
scatter(t_GEO,ic2_GEO,'filled','Marker','v');
ylim([0,42]);
ylabel('Number of iterations','FontName','Times New Roman');
xlabel('Ordinal number of step','FontName','Times New Roman');
legend("$dt$","$I_{c1} $","$I_{c2} $", 'Interpreter', 'latex');