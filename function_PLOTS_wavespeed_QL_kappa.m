function function_PLOTS_wavespeed_QL_kappa

%% Plots for wavespeed and QL as a function of kappa

%% Define phi
phi=1;

%% Data from repeated simulations using code used to generate Figure 1.

kappa_vector = [0.0100
    0.1000
    0.2000
    0.3000
    0.4000
    0.5000
    0.6000
    0.7000
    0.8000
    0.9000
    1.0000
    1.1000
    1.2000
    1.3000
    1.4000
    1.5000
    1.6000
    1.7000
    1.8000
    1.9000
    2.0000];

wavespeed_vector =[ -0.5220
   -0.4720
   -0.4170
   -0.3630
   -0.3090
   -0.2560
   -0.2040
   -0.1520
   -0.1010
   -0.0500
         0
    0.0500
    0.0990
    0.1480
    0.1970
    0.2450
    0.2940
    0.3420
    0.3890
    0.4370
    0.4840];


 qL_vector = [1.8720
    1.7470
    1.6210
    1.5090
    1.4110
    1.3230
    1.2440
    1.1730
    1.1100
    1.0520
    1.0000
    0.9520
    0.9080
    0.8680
    0.8310
    0.7970
    0.7650
    0.7360
    0.7090
    0.6830
    0.6610];
   
filepath_save_figs = pwd;

%% Leading order perturbation solution

len_param_ind_to_plot = length(kappa_vector); %number of kappa values included

wavespeed_vector_leadingperturbation = zeros(len_param_ind_to_plot,1);
qL_vector_leadingperturbation = zeros(len_param_ind_to_plot,1);

% determine the wavespeed and qL using Eqs. (12) and (10), respectively
for ii=1:len_param_ind_to_plot
    kappa_ii=kappa_vector(ii);
    
    c= function_wavespeed_leadingorderperturbation(kappa_ii,phi);
    qL =  function_qL_exact(kappa_ii,phi,c);
    wavespeed_vector_leadingperturbation(ii) = c;
    qL_vector_leadingperturbation(ii) = qL;
end

%% Plot the solutions from the PDE vs the leading order perturbation solution on a combined axis.

figure
yyaxis left
plot(kappa_vector, wavespeed_vector) % PDE solutions
xlabel('Kappa')
ylabel('Wavespeed')
hold on
plot(kappa_vector, wavespeed_vector_leadingperturbation) %Leading order perturbation solution

hold on
yyaxis right
plot(kappa_vector, qL_vector) % PDE solutions
ylabel('qL')
hold on
plot(kappa_vector,qL_vector_leadingperturbation) %Leading order perturbation solution
xlim([0,2])
xticks([0,1,2])
ylim([0,2.5])
yticks([0,1,2])

title('Wavespeed and QL as a function of kappa')
legend('PDE Solution','Leading perturbation solution')

shg

print(gcf,'-depsc2',[filepath_save_figs '\' 'c_QL_vary_kappa_comparison' '.eps'])
saveas(gcf,[filepath_save_figs '\' 'c_QL_vary_kappa_comparison' '.fig'])
saveas(gcf,[filepath_save_figs '\' 'c_QL_vary_kappa_comparison' '.jpg'])



end