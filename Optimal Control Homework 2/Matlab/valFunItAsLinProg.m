function [V,U] = valFunItAsLinProg(F,F0,alpha)
    % VALFUNITASLINPROG Value function iteration. 
    %   [V,U] = VALFUNITASLINPROG(F,F0,alpha) returns solution and
    %   optimal feedback for a discrete-time, infinite horizon DP via a
    %   linear program interpretation. F is the transition matrix
    %   (considering state and control), F0 the cost matrix and alpha the
    %   dicounted cost factor.
    
    [m, n] = size(F);
    
    % construct inequality constraint matricies
    A = repmat(eye(n),m,1) - alpha*onehot(reshape(F',numel(F),1),n);
    b = reshape(F0',numel(F0),1);
    c = -ones(n,1);
    
    % solve problem
    V = linprog(c,A,b);
    
    % Get the optimal feedback
    [~,U] = min(abs(reshape(A*V - b,n,m)),[],2);
    U = U-1;
    
end