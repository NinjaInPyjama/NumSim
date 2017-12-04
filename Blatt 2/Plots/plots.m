close all
clear
clc

DATA = xlsread("tabelle.xlsx");
c = DATA(:,1);
comp = DATA(:,2);
comm = DATA(:,3);
dt = DATA(:,4);
res = DATA(:,5);
runT = DATA(1:7,6);

% figure(1)
% hold on
% scatter(c,comp,'filled');
% scatter(c,comm,'filled');
% hold off
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% grid on
% title('Double logarithmic plot for computation and communication time on 420x420 grid');
% legend('Computation Time', 'Communication Time');
% xlabel('Amount processes');
% ylabel('Time [ms]');
% 
% figure(2)
% hold on
% scatter(c,dt,'filled');
% scatter(c,res,'filled');
% hold off
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% grid on
% title('Double logarithmic plot of time step and global residual calculation on 420x420 grid');
% legend('Time Step Computation', 'Global Residual Computation');
% xlabel('Amount processes');
% ylabel('Time [ms]');

figure(3)
hold on
scatter(c(1:length(runT)),runT,'filled');
hold off
set(gca,'xscale','log')
set(gca,'yscale','log')
grid on
title('Double logarithmic plot of run time on 120x120 grid');
xlabel('Amount processes');
ylabel('Time [s]');

eff = runT(1)./(c(1:length(runT)).*runT);
figure(4)
hold on
scatter(c(1:length(eff)),eff,'filled');
hold off
set(gca,'xscale','log')
set(gca,'yscale','log')
grid on
title('Double logarithmic plot of parallel efficiency on 120x120 grid');
xlabel('Amount processes');
ylabel('Efficiency');
