function [L] = function_ctm_boundary_position(L_old,q,dt,dz,nodesz)

    L = L_old - dt*1/(q(nodesz)^3*L_old*dz)*(  q(nodesz) - q(nodesz-1)) ;

end