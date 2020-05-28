function [g] = function_discretised_func_q(q,q_old,nodesz,dt,dz,L_old,L,kappa,phi)

z=0:dz:1;

r = dt/(dz^2*L_old^2);
dL_dt = (L-L_old)/dt;
g = zeros(nodesz,1);

%% First node
g(1) = q(2)-q(1); 

%% Interior nodes
     for i = 2:nodesz-1
        
        C = z(i)*dt/(2*L_old*dz);
        
        g(i) = - q(i) + q_old(i) ...
            -r*(1/q(i-1) - 2*1/q(i) + 1/q(i+1))...
            + C*dL_dt*(q(i+1)-q(i-1))...
            ...%+ dt*q_old(i)*(1-q_old(i));
            + dt*q(i)*(1-q(i));
    end


%% Final node
i = nodesz;
g(i) = (kappa- 1/q(i)) + (1/(L_old*dz))*(phi/(q(i)^3))*(q(i) - q(i-1));



end