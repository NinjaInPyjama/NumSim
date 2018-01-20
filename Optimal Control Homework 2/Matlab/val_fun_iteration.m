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
    %   realtiv error tolerance and V0 sets the inital state for the fix
    %   point iteration.
    
    if(nargin == 5) 
        V = V0;
    else
        V = zeros(1,8);
        epsi = 1e-16;
    end
    
    U = zeros(1,8);
    
    err = 2*epsi;
    count = 0;
    
    while(err > epsi)
        count = count+1;
        [V_new, U] = min(F0 + alpha*V(F));
        err = norm(V - V_new);
        V = V_new;
    end
    
    U = U-1;
end