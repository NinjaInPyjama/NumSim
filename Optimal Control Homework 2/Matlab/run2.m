%% Setup

clear 
clc

%% Variables

x0 = [0.6; -0.7]; % initial condition in e) and g)
%x0 = [1; -0.9]; % initial condition in f)
A = [1 3; -0.5 1];
P = [4.2 7; 7 36.1];
B = [0; 1];
K = [-0.3, 1.4];
c = min(eig(P))/norm(K).^2;

%% Matricies for QP

H = eye(11);
H(7:8,7:8) = P;

A_eq = -eye(8,11);
A_eq(1:2,1:2) = eye(2);
A_eq = A_eq + [kron(diag(ones(1,3),-1),A),zeros(8,3)] + ...
        [zeros(8),[zeros(2,3);kron(eye(3),B)]];
    
B_eq = [x0; zeros(6,1)];

A_ineq = [zeros(6,8),[eye(3);-eye(3)]];
B_ineq = ones(1,6);

T = zeros(11);
T(7:8,7:8) = P;

d = c;

%% MPC

X_MPC = zeros(2,31);
X_MPC(:,1) = x0;
U_MPC = zeros(1,30);

for k=1:30
    x_start = X_MPC(:,k);
    B_eq = [x_start; zeros(6,1)];
    x_pred = fmincon(@(z) z'*H*z, [x_start; zeros(9,1)], A_ineq,B_ineq,A_eq,B_eq,[],[],@(z) quadCon(T,c,z));
    X_MPC(:,k+1) = x_pred(3:4);
    U_MPC(k) = x_pred(9);
end

%% Plot MPC

f = @(X) diag(X'*P*X);

x1 = linspace(-1.5,1.5,100);
x2 = linspace(-0.5,0.5,100);

[X1,X2] = meshgrid(x1,x2);
Z = reshape(f([X1(:),X2(:)]'),100,100);

figure
set(gca,'fontsize',16);
hold on
plot(X_MPC(1,:),X_MPC(2,:),'color',[0 114 189]/255,'LineWidth',2);
contour(X1,X2,Z,[c,c],'color',[217 83 25]/255,'LineWidth',2);
scatter(X_MPC(1,:),X_MPC(2,:),100,'markeredgecolor',[0 114 189]/255,'LineWidth',1);
hold off
xlim([-1.8, 1.5])
xlabel('$x_1$','Interpreter','latex','FontSize',24);
ylabel('$x_2$','Interpreter','latex','FontSize',24);
L1 = legend('$x_k$ trajectory', 'terminal region');
set(L1,'Interpreter','latex');
set(L1,'FontSize',24);
set(L1,'Location','northeast');
grid on

figure
set(gca,'fontsize',16);
hold on
plot(0:length(U_MPC)-1,U_MPC,'color',[0 114 189]/255,'LineWidth',2);
plot([0, length(U_MPC)-1],[1,1],'--','color',[217 83 25]/255,'LineWidth',2);
scatter(0:length(U_MPC)-1,U_MPC,'markeredgecolor',[0 114 189]/255,'LineWidth',1);
plot([0, length(U_MPC)-1],[-1,-1],'--','color',[217 83 25]/255,'LineWidth',2);
hold off
ylim([-1.5 1.5]);
xlim([0 length(U_MPC)-1])
xlabel('$k$','Interpreter','latex','FontSize',24);
ylabel('$u_k$','Interpreter','latex','FontSize',24);
L1 = legend('$x_k$ trajectory', 'terminal region');
set(L1,'Interpreter','latex');
set(L1,'FontSize',24);
set(L1,'Location','northeast');
grid on

%% LQ

X_LQ = zeros(2,31);
X_LQ(:,1) = x0;
U_LQ = zeros(1,30);

[K,~,~] = dlqr(A,B,eye(2),1,0);

for k=1:30
    U_LQ(k) = -K*X_LQ(:,k);
    X_LQ(:,k+1) = (A - B*K)*X_LQ(:,k);
end

%% Plot LQ

figure
set(gca,'fontsize',16);
hold on
plot(X_LQ(1,:),X_LQ(2,:),'color',[0 114 189]/255,'LineWidth',2);
contour(X1,X2,Z,[c,c],'color',[217 83 25]/255,'LineWidth',2);
scatter(X_LQ(1,:),X_LQ(2,:),100,'markeredgecolor',[0 114 189]/255,'LineWidth',1);
hold off
xlim([-1.8, 1.5])
xlabel('$x_1$','Interpreter','latex','FontSize',24);
ylabel('$x_2$','Interpreter','latex','FontSize',24);
L1 = legend('$x_k$ trajectory', 'terminal region');
set(L1,'Interpreter','latex');
set(L1,'FontSize',24);
set(L1,'Location','northeast');
grid on

figure
set(gca,'fontsize',16);
hold on
plot(0:length(U_LQ)-1,U_LQ,'color',[0 114 189]/255,'LineWidth',2);
scatter(0:length(U_LQ)-1,U_LQ,100,'markeredgecolor',[0 114 189]/255,'LineWidth',1);
hold off
ylim([-1.5 1.5]);
xlim([0 length(U_LQ)-1])
xlabel('$k$','Interpreter','latex','FontSize',24);
ylabel('$u_k$','Interpreter','latex','FontSize',24);
grid on
 