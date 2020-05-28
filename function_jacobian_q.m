function [L_diag, D_diag, U_diag] = function_jacobian_q(q,nodesz,dt,dz,L_old,L,kappa, phi)

z=0:dz:1;

L_diag = zeros(nodesz,1);
D_diag = zeros(nodesz,1);
U_diag = zeros(nodesz,1);

r = dt/(dz^2*L_old^2);

dL_dt = (L-L_old)/dt;


%% First node
D_diag(1) = -1;
U_diag(1) = 1;

%% Internal nodes

for i = 2:nodesz-1
    C = z(i)*dt/(2*L_old*dz);
    D_diag(i) = -1 - 2*r*1/q(i)^2 + dt*(1-2*q(i));
    L_diag(i) = r*1/q(i-1)^2 - C*dL_dt;
    U_diag(i) = r*1/q(i+1)^2 + C*dL_dt;
end


%% Final node

i = nodesz;


D_diag(i) =   1/(q(i)^2) + (1/(L_old*dz))*(phi/(q(i)^3)) + (1/(L_old*dz))*phi*(-3/(q(i)^4))*(q(i) - q(i-1));
L_diag(i) = (1/(L_old*dz))*(phi/(q(i)^3))*(-1);



end