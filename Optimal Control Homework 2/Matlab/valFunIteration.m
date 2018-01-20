function [V,U,count] = valFunIteration(F,F0,alpha,epsi,V0)
    % VALFUNITERATION Value function iteration. 
    %   [V,U,count] = VALFUNITERATION(F,F0,alpha) returns solution and
    %   optimal feedback for a discrete-time (and iteration count of the
    %   fix point iteration, infinite horizon DP via value function
    %   iteration. F is the transition matrix (considering state and
    %   control), F0 the cost matrix and alpha dicounted cost factor. (Sets
    %   epsi = 1e-16 and V = 0).
    %
    %   [V,U,count] = VALFUNITERATION(F,F0,alpha,epsi,V0) sets the
    %   error tolerance and V0 sets the inital state for the fix point
    %   iteration.
    
    n = size(F,2);
    
    if(nargin == 5) 
        V = V0;
    else
        V = zeros(1,n); % inital value
        epsi = 1e-16; % error tolerance
    end
    
    err = 2*epsi; %% current error (pre-set)
    count = 0; % amount of iterations for the fix point iteration
    
    % Fix point iteration
    while(err > epsi)
        count = count+1; % increase counter
        [V_new, U] = min(F0 + alpha*V(F)); % calculate V_{k+1} = TV_k
        err = norm(V - V_new); % calculate new error
        V = V_new; % Update V
    end
    
    V = V';
    U = U'-1;
end