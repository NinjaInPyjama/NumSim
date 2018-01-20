clear 
clc

x0 = [0.6; -0.7]; %[1; -0.9]; %
A = [1 3; -0.5 1];
P = [4.2 7; 7 36.1];
B = [0; 1];
K = [-0.3, 1.4];
c = min(eig(P))/norm(K).^2;


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

X = zeros(2,31);
X(:,1) = x0;
U = zeros(1,3);

% for k=1:30
%     x_start = X(:,k);
%     B_eq = [x_start; zeros(6,1)];
%     x_pred = fmincon(@(z) z'*H*z, [x_start; zeros(9,1)], A_ineq,B_ineq,A_eq,B_eq,[],[],@(z) quadCon(T,c,z));
%     X(:,k+1) = x_pred(3:4);
%     U(k) = x_pred(9);
% end

[K,~,~] = dlqr(A,B,eye(2),1,0);

for k=1:30
    U(k) = -K*X(:,k);
    X(:,k+1) = (A - B*K)*X(:,k);
end

figure
plot(0:30,X(1,:),0:30,X(2,:),0:29,U,'LineWidth',2)

