function function_TWAVE_pdesolution(simulation_id,L,timestop,time_record_vector_ctm,kappa,phi)

%% PDE solution

%Discretisations and Newton-Raphson parameters

dz=0.00001; %for accuracy
%dz=0.0001; %for speed, solution is also accurate
nodesz=round(1/dz) +1;
dt=1e-2; %timestep

err_tol = 1e-8;
max_iters=1000;

%% q initial condition

function_q_init_condition_9 = @(x) 0.*x  + 1;

x_disc = (0:dz:1)*L;

q0 = function_q_init_condition_9(x_disc)';


%% create a folder to save the results

file_save_name = ['Results_Sim' simulation_id '_Continuum'];

folder_name = [simulation_id '_Continuum'];
if ~exist([folder_name'], 'dir')
    mkdir([folder_name]);
end

%% convert time_record_vector_ctm into timesteps

t_steprec=0;

time_record_vector_dt = round(time_record_vector_ctm/dt,0);


%% storage and save initial values

loop_count_stored = 1;
max_recorded_data_points = length(time_record_vector_dt); %value to initialise for below

q_hist = zeros(nodesz,max_recorded_data_points);
q_hist(:,loop_count_stored) = q0;
q=q0;

t_hist  = zeros(1,max_recorded_data_points);
t_hist(1)=0;

L_hist   = zeros(1,max_recorded_data_points);
L_hist(1)=L;

%% Run the Newton-Raphson method to numerically solve the PDE

t=0; %initial time

q_update = q;
L_update = L;

while t < timestop
    
    t = t + dt; %update time
    
    q_old = q;
    L_old = L;
    
    res = err_tol + 1; %intialise the residual to be updated inside the Newton-Raphson iteration for this timestep
    
    w=1; %initialise the iteration counter for this timestep
    
    %Newton-Raphson iterations
    while w <= max_iters && res > err_tol
        
        if L_update < 1e-2 %for negative wavespeeds stop solution close to boundary.
            w = max_iters + 1;
        else
            
            
            q = q_update;
            L = L_update;
            
            
            
            %% Solve for q
            [L_diag,D_diag,U_diag] = function_jacobian_q(q,nodesz,dt,dz,L_old,L,kappa,phi);
            rhs = -function_discretised_func_q(q,q_old,nodesz,dt,dz,L_old,L,kappa,phi);
            dq = function_tridia(nodesz,L_diag,D_diag,U_diag,rhs); % solve solution using thomas algorthim
            q_update = q + dq;
            
            
            %% Solve for L
            [L_update] = function_ctm_boundary_position(L_old,q,dt,dz,nodesz);
            
            
            %% Calculate the residual
            res = norm(dq,inf); % update residual
            
            w=w+1;
            
        end
        
    end
    
    % If exceed the maximum number of iterations output a message.
    if w >= max_iters
        fprintf('\n MAX ITERS REACHED! :t = %f; iters = %f\n', t,w);
    end
    
    
    %% Store data if near a timestep where data should be stored.
    
    t_steprec = t_steprec + 1;
    
    if sum(time_record_vector_dt == t_steprec) > 0
        
        loop_count_stored = loop_count_stored + 1;
        q_hist(:,loop_count_stored) = q_update;
        t_hist(loop_count_stored) = t;
        L_hist(loop_count_stored) = L_update;
        
        fprintf('\nt = %f; iters = %f\n', t,w);
        
    end
    
    %Error message if density blows up (shouldnt happen and doesnt for parameters considered)
    if sum(isnan(q_update)) > 1
        fprintf('\n q is NAN t = %f; iters = %f\n', t,w);
        break
    end
    
    
    
end


%% Save the full .mat file.

save([pwd '\' folder_name '\' file_save_name],'-v7.3');


end