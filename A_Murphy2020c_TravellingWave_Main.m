%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2020--Murphy et al--Travelling waves in a free boundary mechanobiological model of an epithelial tissue

% Key algorithms used to generate paper Figures 1b,d and Figure 2:

% Figure 1b,d for kappa=2,phi=1.
% Figure 1b - pde solution to travelling wave with postive speed.
% Figure 1d - phase plane for travelling wave with positive speed.

% Figure 2 for 0 < kappa < 2, phi=1.
% Figure 2 - Wavespeed as a function of kappa, boundary density as a function of kappa, leading order perturbation solution vs pde solution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Reset matlab
clear
clc

%% Figure 1b,d and leading order perturbation solution vs pde solution.

for parameter_id =3001 % choose the simulation id
    
    for ctm_plot=1:2 %generate the pde solution and plot
        
        if parameter_id ==3001
             
            simulation_id = 'GitHub_TravelWave_3001';
            L=10; %initial tissue length
            timestop=1000; 
            time_record_vector_ctm=0:5:timestop; %when to record numerical solution
            kappa=2;
            phi=1;
            
        end
        
        if   ctm_plot == 1
            %% PDE solution
            
            % Eqs. (1)-(4) in the paper can be solved numerically can be solved numerically by using a boundary fixing transformation \cite{Kutluay1997},
            % discretising the subsequent equations on a uniform mesh using a central difference approximation. 
            % The resulting system of ordinary differential equations are solved using an implicit Euler approximation, 
            % leading to a system of nonlinear algebraic equations that are solved using Newton-Raphson iteration.
            
            function_TWAVE_pdesolution(simulation_id,L,timestop,time_record_vector_ctm,kappa,phi)
            
        elseif ctm_plot == 2
            %% Plots
            % Density snapshots
            % Phase plane
            % Leading order perturbation solution v PDE solution
            function_TWAVE_PLOTS(simulation_id)
            
            
        end

    end
end  


%% Figure 2 - wavespeed as a function of kappa, Q_L as a function of kappa
% Results for 0 < kappa < 2 have been obtained by looping over kappa in the code used to generate Figure 1.
function_PLOTS_wavespeed_QL_kappa
    
