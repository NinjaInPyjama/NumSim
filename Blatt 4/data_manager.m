%%
%   _____  ___   ____  ____  ___      ___   ________  __     ___      ___ 
%  (\"   \|"  \ ("  _||_ " ||"  \    /"  | /"       )|" \   |"  \    /"  |
%  |.\\   \    ||   (  ) : | \   \  //   |(:   \___/ ||  |   \   \  //   |
%  |: \.   \\  |(:  |  | . ) /\\  \/.    | \___  \   |:  |   /\\  \/.    |
%  |.  \    \. | \\ \__/ // |: \.        |  __/  \\  |.  |  |: \.        |
%  |    \    \ | /\\ __ //\ |.  \    /:  | /" \   :) /\  |\ |.  \    /:  |
%   \___|\____\)(__________)|___|\__/|___|(_______/ (__\_|_)|___|\__/|___|

% AUTHORS: FEINLER, Mathias
%          GRAMLING, Silvia
%          SCHMIDGALL, Markus

%% GENERAL SETTINGS
clear
clc

%% LOADING DATA (MONTE CARLO)

% amount of simulations ran
num_simulations_mc = 2016;

% reynolds Numbers
re_mc = zeros(1,num_simulations_mc);

% load first simulation data
buffer = load("MonteCarlo030118/output_0.txt");

% store first reynolds number
re_mc(1) = buffer(1,1);

% store time steps
buffer(1,1) = .0;
t_mc = buffer(:,1);
num_t_steps_mc = size(t_mc,1);

% alocate and store simulation data
data_mc = zeros(size(buffer,1),3*num_simulations_mc);
data_mc(:,1:3) = buffer(:,2:4);

for i = 1:num_simulations_mc-1
    % load simulation data
    buffer = load(['MonteCarlo030118/output_', num2str(i), '.txt']);
    
    % store first reynolds number
    re_mc(i+1) = buffer(1,1);

    % store simulation data
    data_mc(:,3*i+(1:3)) = buffer(:,2:4);
end


%% PROCESS DATA (MONTE CARLO)

% alocate expectation and variance of velocities
E_mc = zeros(num_t_steps_mc,3);
V_mc = zeros(num_t_steps_mc,3);

for i=0:2
    U = data_mc(:,i+(1:3:3*num_simulations_mc));
    E_mc(:,i+1) = sum(U,2)/num_simulations_mc;
    V_mc(:,i+1) = sum((U - E_mc(:,i+1)).^2,2)/(num_simulations_mc-1);
end

%% CONVERGENCE STUDY (MONTE CARLO)

INDEX = (1:num_simulations_mc)';%randperm(num_simulations_mc)';
INDEX = reshape([3*(INDEX-1)+1, 3*(INDEX-1)+2, 3*INDEX]',3*num_simulations_mc,1); 
DATA = data_mc(:,INDEX);

E_mc_conv = zeros(num_simulations_mc-1,3);
V_mc_conv = zeros(num_simulations_mc-1,3);

for k=2:num_simulations_mc
    for i=0:2
        U = DATA(end,i+(1:3:3*k));
        E_mc_conv(k-1,i+1) = sum(U,2)/k;
        V_mc_conv(k-1,i+1) = sum((U - E_mc_conv(k-1,i+1)).^2,2)/(k-1);
    end
end

E_mc_conv = abs(E_mc_conv - E_mc_conv(end,:));
V_mc_conv = abs(V_mc_conv - V_mc_conv(end,:));

f = figure;
set(gca,'FontSize',10);
p = uipanel('Parent',f,'BorderType','none'); 
p.Title = 'Convergence Study for Monte Carlo Method'; 
p.TitlePosition = 'centertop'; 
p.FontSize = 20;
p.FontWeight = 'bold';

subplot(3,2,1,'Parent',p) 
loglog(2:num_simulations_mc, E_mc_conv(:,1),[500,1000], E_mc_conv([499,999],1),'linewidth',2)
grid on
title('Error of Expectation of $u_{120,5}$','Interpreter','latex','FontSize',14);
xlabel('$Re$ numbers considered','Interpreter','latex');
ylabel('error','Interpreter','latex');

subplot(3,2,2,'Parent',p) 
loglog(2:num_simulations_mc, V_mc_conv(:,1),[500,1000], V_mc_conv([499,999],1),'linewidth',2)
grid on
title('Error of Variation of $u_{120,5}$','Interpreter','latex','FontSize',14);
xlabel('$Re$ numbers considered','Interpreter','latex');
ylabel('error','Interpreter','latex');

subplot(3,2,3,'Parent',p) 
loglog(2:num_simulations_mc, E_mc_conv(:,2),[500,1000], E_mc_conv([499,999],2),'linewidth',2)
grid on
title('Error of Expectation of $u_{64,64}$','Interpreter','latex','FontSize',14);
xlabel('$Re$ numbers considered','Interpreter','latex');
ylabel('error','Interpreter','latex');

subplot(3,2,4,'Parent',p) 
loglog(2:num_simulations_mc, V_mc_conv(:,2),[500,1000], V_mc_conv([499,999],2),'linewidth',2)
grid on
title('Error of Variation of $u_{64,64}$','Interpreter','latex','FontSize',14);
xlabel('$Re$ numbers considered','Interpreter','latex');
ylabel('error','Interpreter','latex');

subplot(3,2,5,'Parent',p) 
loglog(2:num_simulations_mc, E_mc_conv(:,3),[500,1000], E_mc_conv([499,999],3),'linewidth',2)
grid on
title('Error of Expectation of $u_{5,120}$','Interpreter','latex','FontSize',14);
xlabel('$Re$  numbers considered','Interpreter','latex');
ylabel('error','Interpreter','latex');

subplot(3,2,6,'Parent',p) 
loglog(2:num_simulations_mc, V_mc_conv(:,3),[500,1000], V_mc_conv([499,999],3),'linewidth',2)
grid on
title('Error of Variation of $u_{5,120}$','Interpreter','latex','FontSize',14);
xlabel('$Re$  numbers considered','Interpreter','latex');
ylabel('error','Interpreter','latex');


%% PLOT DATA (MONTE CARLO)

f = figure;
set(gca,'FontSize',10);
p = uipanel('Parent',f,'BorderType','none'); 
p.Title = 'Simulation results for Monte Carlo Method'; 
p.TitlePosition = 'centertop'; 
p.FontSize = 20;
p.FontWeight = 'bold';

subplot(3,2,1,'Parent',p) 
plot(t_mc,E_mc(:,1),'linewidth',3)
grid on
title('Expectation of $u_{120,5}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');

subplot(3,2,2,'Parent',p)
plot(t_mc ,V_mc(:,1),'linewidth',3)
grid on
title('Variance of $u_{120,5}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');

subplot(3,2,3,'Parent',p) 
plot(t_mc,E_mc(:,2),'linewidth',3)
grid on
title('Expectation of $u_{64,64}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');

subplot(3,2,4,'Parent',p) 
plot(t_mc,V_mc(:,2),'linewidth',3)
grid on
title('Variance of $u_{64,64}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');

subplot(3,2,5,'Parent',p) 
plot(t_mc,E_mc(:,3),'linewidth',3)
grid on
title('Expectation of $u_{5,120}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');

subplot(3,2,6,'Parent',p) 
plot(t_mc,V_mc(:,3),'linewidth',3)
grid on
title('Variance of $u_{5,120}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');

mean_re = sum(re_mc)/num_simulations_mc;
var_re = sqrt(sum((re_mc - mean_re).^2)/(num_simulations_mc - 1));
my_re = 1000:2000;
figure
hold on
plot(my_re,normpdf(my_re,1500,1000/6),'--','linewidth',2)
plot(my_re,normpdf(my_re,mean_re,var_re),'linewidth',3)
hold off
set(gca,'FontSize',14);
grid on
title('Approximated distribution of Reynolds numbers for Monte Carlo Method','Interpreter','latex','FontSize',24);
xlabel('$Re$','Interpreter','latex');
L = legend('Exact distribution','Approximated distribution');
set(L,'Interpreter','latex');
set(L,'FontSize',18);

% figure
% plot(t_mc,E_mc,'linewidth',3)
% set(gca,'FontSize',14);
% grid on
% title('Expectation and of velocities over time','Interpreter','latex','FontSize',24);
% xlabel('$t$ in $s$','Interpreter','latex');
% ylabel('$y$','Interpreter','latex');
% L = legend('$E(u_{120,5})$', '$E(u_{64,64})$', '$E(u_{5,120})$');
% set(L,'Interpreter','latex');
% set(L,'FontSize',18);
% 
% figure
% plot(t_mc,V_mc,'linewidth',3)
% set(gca,'FontSize',14);
% grid on
% title('Variance and of velocities over time','Interpreter','latex','FontSize',24);
% xlabel('$t$ in $s$','Interpreter','latex');
% ylabel('$y$','Interpreter','latex');
% L = legend('$V(u_{120,5})$', '$V(u_{64,64})$', '$V(u_{5,120})$');
% set(L,'Interpreter','latex');
% set(L,'FontSize',18);


%% LOADING DATA (TRAPEZ)

% amount simulations ran
num_simulations_tr = 200;

% reynolds Numbers
re_tr = zeros(1,num_simulations_tr);

% load first simulation data
buffer = load("output_trapez/trapez_0.txt");

% store first reynolds number
re_tr(1) = buffer(1,1);

% store time steps
buffer(1,1) = .0;
t_tr = buffer(:,1);
t_tr(end) = 15.0;
num_t_steps_tr = size(t_tr,1);

% alocate and store simulation data
data_tr = zeros(size(buffer,1),3*num_simulations_tr);
data_tr(:,1:3) = buffer(:,2:4);

for i = 1:num_simulations_tr-1
    % load simulation data
    buffer = load(['Trapezregel140118/output_', num2str(i), '.txt']);
    
    % store first reynolds number
    re_tr(i+1) = buffer(1,1);

    % store simulation data
    data_tr(:,3*i+(1:3)) = buffer(:,2:4);
end


%% PROCESS DATA (TRAPEZ)

% data matices for different amount of nodes
data_tr_10 = data_tr(:,reshape(repmat(3*(1:20:num_simulations_tr)-3,3,1)+(1:3)',1,3*ceil(num_simulations_tr/20)));
data_tr_25 = data_tr(:,reshape(repmat(3*(1:8:num_simulations_tr)-3,3,1)+(1:3)',1,3*ceil(num_simulations_tr/8)));
data_tr_50 = data_tr(:,reshape(repmat(3*(1:4:num_simulations_tr)-3,3,1)+(1:3)',1,3*ceil(num_simulations_tr/4)));
data_tr_100 = data_tr(:,reshape(repmat(3*(1:2:num_simulations_tr)-3,3,1)+(1:3)',1,3*ceil(num_simulations_tr/2)));
data_tr_200 = data_tr;

% alocate expectation and variance of velocities (50, 100 and 200 nodes)
E_tr_10 = zeros(num_t_steps_tr,3);
V_tr_10 = zeros(num_t_steps_tr,3);
E_tr_25 = zeros(num_t_steps_tr,3);
V_tr_25 = zeros(num_t_steps_tr,3);
E_tr_50 = zeros(num_t_steps_tr,3);
V_tr_50 = zeros(num_t_steps_tr,3);
E_tr_100 = zeros(num_t_steps_tr,3);
V_tr_100 = zeros(num_t_steps_tr,3);
E_tr_200 = zeros(num_t_steps_tr,3);
V_tr_200 = zeros(num_t_steps_tr,3);

for i=0:2
    U = data_tr_10(:,i+(1:3:size(data_tr_10,2)));
    prop = normpdf(repmat(re_tr(1:20:num_simulations_tr),num_t_steps_tr,1),1500,1000/6);
    E_tr_10(:,i+1) = 100*sum(U.*prop./repmat([.5,ones(1,size(data_tr_10,2)/3-2), .5],num_t_steps_tr,1),2);
    V_tr_10(:,i+1) = 100*sum((U - E_tr_10(:,i+1)).^2.*prop./repmat([.5,ones(1,size(data_tr_10,2)/3-2), .5],num_t_steps_tr,1),2);
    
    U = data_tr_25(:,i+(1:3:size(data_tr_25,2)));
    prop = normpdf(repmat(re_tr(1:8:num_simulations_tr),num_t_steps_tr,1),1500,1000/6);
    E_tr_25(:,i+1) = 40*sum(U.*prop./repmat([.5,ones(1,size(data_tr_25,2)/3-2), .5],num_t_steps_tr,1),2);
    V_tr_25(:,i+1) = 40*sum((U - E_tr_25(:,i+1)).^2.*prop./repmat([.5,ones(1,size(data_tr_25,2)/3-2), .5],num_t_steps_tr,1),2);
    
    U = data_tr_50(:,i+(1:3:size(data_tr_50,2)));
    prop = normpdf(repmat(re_tr(1:4:num_simulations_tr),num_t_steps_tr,1),1500,1000/6);
    E_tr_50(:,i+1) = 20*sum(U.*prop./repmat([.5,ones(1,size(data_tr_50,2)/3-2), .5],num_t_steps_tr,1),2);
    V_tr_50(:,i+1) = 20*sum((U - E_tr_50(:,i+1)).^2.*prop./repmat([.5,ones(1,size(data_tr_50,2)/3-2), .5],num_t_steps_tr,1),2);
    
    U = data_tr_100(:,i+(1:3:size(data_tr_100,2)));
    prop = normpdf(repmat(re_tr(1:2:num_simulations_tr),num_t_steps_tr,1),1500,1000/6);
    E_tr_100(:,i+1) = 10*sum(U.*prop./repmat([.5,ones(1,size(data_tr_100,2)/3-2), .5],num_t_steps_tr,1),2);
    V_tr_100(:,i+1) = 10*sum((U - E_tr_100(:,i+1)).^2.*prop./repmat([.5,ones(1,size(data_tr_100,2)/3-2), .5],num_t_steps_tr,1),2);
    
    U = data_tr_200(:,i+(1:3:size(data_tr_200,2)));
    prop = normpdf(repmat(re_tr,num_t_steps_tr,1),1500,1000/6);
    E_tr_200(:,i+1) = 5*sum(U.*prop./repmat([.5,ones(1,size(data_tr_200,2)/3-2), .5],num_t_steps_tr,1),2);
    V_tr_200(:,i+1) = 5*sum((U - E_tr_200(:,i+1)).^2.*prop./repmat([.5,ones(1,size(data_tr_200,2)/3-2), .5],num_t_steps_tr,1),2);
end

%% CONVERGENCE STUDY (TRAPEZ)

E_tr_conv = abs([E_tr_10(end,:);E_tr_25(end,:);E_tr_50(end,:);E_tr_100(end,:)] - E_tr_200(end,:));
V_tr_conv = abs([V_tr_10(end,:);V_tr_25(end,:);V_tr_50(end,:);V_tr_100(end,:)] - V_tr_200(end,:));

f = figure;
set(gca,'FontSize',10);
p = uipanel('Parent',f,'BorderType','none'); 
p.Title = 'Convergence Study for Trapez Method'; 
p.TitlePosition = 'centertop'; 
p.FontSize = 20;
p.FontWeight = 'bold';

subplot(3,2,1,'Parent',p) 
loglog([10,25,50,100], E_tr_conv(:,1),'linewidth',2)
grid on
title('Error of Expectation of $u_{120,5}$','Interpreter','latex','FontSize',14);
xlabel('$Re$ numbers considered','Interpreter','latex');
ylabel('error','Interpreter','latex');

subplot(3,2,2,'Parent',p) 
loglog([10,25,50,100], V_tr_conv(:,1),'linewidth',2)
grid on
title('Error of Variation of $u_{120,5}$','Interpreter','latex','FontSize',14);
xlabel('$Re$ numbers considered','Interpreter','latex');
ylabel('error','Interpreter','latex');

subplot(3,2,3,'Parent',p) 
loglog([10,25,50,100], E_tr_conv(:,2),'linewidth',2)
grid on
title('Error of Expectation of $u_{64,64}$','Interpreter','latex','FontSize',14);
xlabel('$Re$ numbers considered','Interpreter','latex');
ylabel('error','Interpreter','latex');

subplot(3,2,4,'Parent',p) 
loglog([10,25,50,100], V_tr_conv(:,2),'linewidth',2)
grid on
title('Error of Variation of $u_{64,64}$','Interpreter','latex','FontSize',14);
xlabel('$Re$ numbers considered','Interpreter','latex');
ylabel('error','Interpreter','latex');

subplot(3,2,5,'Parent',p) 
loglog([10,25,50,100], E_tr_conv(:,3),'linewidth',2)
grid on
title('Error of Expectation of $u_{5,120}$','Interpreter','latex','FontSize',14);
xlabel('$Re$  numbers considered','Interpreter','latex');
ylabel('error','Interpreter','latex');

subplot(3,2,6,'Parent',p) 
loglog([10,25,50,100], V_tr_conv(:,3),'linewidth',2)
grid on
title('Error of Variation of $u_{5,120}$','Interpreter','latex','FontSize',14);
xlabel('$Re$  numbers considered','Interpreter','latex');
ylabel('error','Interpreter','latex');

%% PLOT DATA (TRAPEZ)

% 50 Nodes
f = figure;
set(gca,'FontSize',10);
p = uipanel('Parent',f,'BorderType','none'); 
p.Title = 'Simulation results for Trapez Rule (50 Nodes)'; 
p.TitlePosition = 'centertop'; 
p.FontSize = 20;
p.FontWeight = 'bold';

subplot(3,2,1,'Parent',p) 
plot(t_tr,E_tr_50(:,1),'linewidth',3)
grid on
title('Expectation of $u_{120,5}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');

subplot(3,2,2,'Parent',p) 
plot(t_tr,V_tr_50(:,1),'linewidth',3)
grid on
title('Variance of $u_{120,5}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');

subplot(3,2,3,'Parent',p) 
plot(t_tr,E_tr_50(:,2),'linewidth',3)
grid on
title('Expectation of $u_{64,64}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');

subplot(3,2,4,'Parent',p) 
plot(t_tr,V_tr_50(:,2),'linewidth',3)
grid on
title('Variance of $u_{64,64}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');

subplot(3,2,5,'Parent',p) 
plot(t_tr,E_tr_50(:,3),'linewidth',3)
grid on
title('Expectation of $u_{5,120}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');

subplot(3,2,6,'Parent',p) 
plot(t_tr,V_tr_50(:,3),'linewidth',3)
grid on
title('Variance of $u_{5,120}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');


% 100 Nodes
f = figure;
set(gca,'FontSize',10);
p = uipanel('Parent',f,'BorderType','none'); 
p.Title = 'Simulation results for Trapez Rule (100 Nodes)'; 
p.TitlePosition = 'centertop'; 
p.FontSize = 20;
p.FontWeight = 'bold';

subplot(3,2,1,'Parent',p) 
plot(t_tr,E_tr_100(:,1),'linewidth',3)
grid on
title('Expectation of $u_{120,5}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');

subplot(3,2,2,'Parent',p) 
plot(t_tr,V_tr_100(:,1),'linewidth',3)
grid on
title('Variance of $u_{120,5}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');

subplot(3,2,3,'Parent',p) 
plot(t_tr,E_tr_100(:,2),'linewidth',3)
grid on
title('Expectation of $u_{64,64}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');

subplot(3,2,4,'Parent',p) 
plot(t_tr,V_tr_100(:,2),'linewidth',3)
grid on
title('Variance of $u_{64,64}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');

subplot(3,2,5,'Parent',p) 
plot(t_tr,E_tr_100(:,3),'linewidth',3)
grid on
title('Expectation of $u_{5,120}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');

subplot(3,2,6,'Parent',p) 
plot(t_tr,V_tr_100(:,3),'linewidth',3)
grid on
title('Variance of $u_{5,120}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');


% 200 Nodes
f = figure;
set(gca,'FontSize',10);
p = uipanel('Parent',f,'BorderType','none'); 
p.Title = 'Simulation results for Trapez Rule (200 Nodes)'; 
p.TitlePosition = 'centertop'; 
p.FontSize = 20;
p.FontWeight = 'bold';

subplot(3,2,1,'Parent',p) 
plot(t_tr,E_tr_200(:,1),'linewidth',3)
grid on
title('Expectation of $u_{120,5}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');

subplot(3,2,2,'Parent',p) 
plot(t_tr,V_tr_200(:,1),'linewidth',3)
grid on
title('Variance of $u_{120,5}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');

subplot(3,2,3,'Parent',p) 
plot(t_tr,E_tr_200(:,2),'linewidth',3)
grid on
title('Expectation of $u_{64,64}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');

subplot(3,2,4,'Parent',p) 
plot(t_tr,V_tr_200(:,2),'linewidth',3)
grid on
title('Variance of $u_{64,64}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');

subplot(3,2,5,'Parent',p) 
plot(t_tr,E_tr_200(:,3),'linewidth',3)
grid on
title('Expectation of $u_{5,120}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');

subplot(3,2,6,'Parent',p) 
plot(t_tr,V_tr_200(:,3),'linewidth',3)
grid on
title('Variance of $u_{5,120}$ over time','Interpreter','latex','FontSize',14);
xlabel('$t$ in $s$','Interpreter','latex');

%% Expectation and variance plots

figure
subplot(2,1,1)
plot(t_tr,E_tr_50(:,2),t_tr,E_tr_100(:,2),t_tr,E_tr_200(:,2),'linewidth',3) %t_tr,E_tr_25(:,2),t_tr,E_tr_10(:,2),
set(gca,'FontSize',14);
grid on
title('Expectation of $u_{64,64}$ over time','Interpreter','latex','FontSize',24);
xlabel('$t$ in $s$','Interpreter','latex');
L = legend('$50$ nodes', '$100$ nodes', '$200$ nodes');
set(L,'Interpreter','latex');
set(L,'FontSize',18);
set(L,'Location','southeast')

subplot(2,1,2)
plot(t_tr,V_tr_50(:,2),t_tr,V_tr_100(:,2),t_tr,V_tr_200(:,2),'linewidth',3) % t_tr,V_tr_25(:,2),t_tr,V_tr_10(:,2),
set(gca,'FontSize',14);
grid on
title('Variance of $u_{64,64}$ over time','Interpreter','latex','FontSize',24);
xlabel('$t$ in $s$','Interpreter','latex');
L = legend('$50$ nodes', '$100$ nodes', '$200$ nodes');
set(L,'Interpreter','latex');
set(L,'FontSize',18);
set(L,'Location','northeast')