function [V,U] = valFunItAsLinProg(F,F0,alpha)
    % VALFUNITASLINPROG Value function iteration. 
    %   [V,U] = VALFUNITASLINPROG(F,F0,alpha) returns solution and
    %   optimal feedback for a discrete-time (and iteration count of the
    %   fix point iteration, infinite horizon DP via  a linear program
    %   interpretation. F is the transition matrix (considering state and
    %   control), F0 the cost matrix and alpha dicounted cost factor. (Sets
    %   epsi = 1e-16 and V = 0).
    %
    %   [V,U] = VALFUNITASLINPROG(F,F0,alpha,epsi,V0) sets the
    %   realtiv error tolerance and V0 sets the inital state for the fix
    %   point iteration.
    
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