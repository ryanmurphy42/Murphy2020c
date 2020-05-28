function function_TWAVE_PLOTS(simulation_id)


%% Load the data

filepath_save_figs = [pwd '\' simulation_id '_Continuum\'];
load([filepath_save_figs 'Results_Sim' simulation_id '_Continuum.mat']);

folder_name = [simulation_id '_Continuum'];
if ~exist([folder_name'], 'dir')
    mkdir([folder_name]);
end

%% Density snapshots

%% Density snapshot - Early times
figure
hold on
plot_times_early=[0,10,20,30,40,50];

for t_loop = 1:length(plot_times_early)
    
    %find the nearest ctm plot time
    [~, first_compare_index]   = min(abs(t_hist(1:end) - plot_times_early(t_loop)));
    
    %overlay the ctm result
    hold on
    plot((0:dz:1)*L_hist(:,first_compare_index) , q_hist(:,first_compare_index),'LineWidth',2)
    
end

legend('t=0', 't=10', 't=20','t=30','t=40','t=50')
xlabel('Position')
ylabel('Density')
title('Density snapshots - early time')
box on
ylim([0,2])
shg

%save figures
print(gcf,'-depsc2',[filepath_save_figs '\' 'Density_evolution_ctm_early' num2str(t_loop) '.eps'])
saveas(gcf,[filepath_save_figs '\' 'Density_evolution_ctm_early' num2str(t_loop) '.fig'])
saveas(gcf,[filepath_save_figs '\' 'Density_evolution_ctm_early' num2str(t_loop) '.jpg'])




%% Density snapshot - Late times

figure
hold on
plot_times_late = [0,200,400,600,800,1000];
for t_loop = 1:length(plot_times_late)
    
    %find the nearest ctm plot time
    [~, first_compare_index]   = min(abs(t_hist(1:end) - plot_times_late(t_loop)));
    
    %overlay the ctm result
    hold on
    plot((0:dz:1)*L_hist(:,first_compare_index) , q_hist(:,first_compare_index),'LineWidth',2)
    
end

%update figure properties

xlabel('Position')
ylabel('Density')
title('Density snapshots - late time')
box on
legend('t=0', 't=200', 't=400','t=600','t=800','t=1000')
ylim([0,2])

print(gcf,'-depsc2',[filepath_save_figs '\' 'Density_evolution_ctm_late' num2str(t_loop) '.eps'])
saveas(gcf,[filepath_save_figs '\' 'Density_evolution_ctm_late' num2str(t_loop) '.fig'])
saveas(gcf,[filepath_save_figs '\' 'Density_evolution_ctm_late' num2str(t_loop) '.jpg'])


%% Plot leading order perturbation solution


figure
hold on

% for positive wavespeeds choose end time
plot_times_perturb = 1000;
% for negative wavespeeds define plot_times_perturb as a time when the travelling wave solution exists (after transient time to form and before boundary effects)

%find the nearest ctm plot time
[~, first_compare_index]   = min(abs(t_hist(1:end) - plot_times_perturb));

%overlay the ctm result
hold on
plot(((0:dz:1)*L_hist(:,first_compare_index) - L_hist(:,first_compare_index)) , q_hist(:,first_compare_index),'LineWidth',2)

c=(L_hist(first_compare_index)-L_hist(first_compare_index-1))/(t_hist(first_compare_index)-t_hist(first_compare_index-1));


% plot the leading order perturbation solution (Eq. 11)
tspan = [0, 10];
y0 = 1/(kappa-c*phi);
if c >0
    [z,Q] = ode45(@(z,Q) Q^2*sqrt(2*(Q-log(Q)-1)), tspan, y0);
elseif c <0
    [z,Q] = ode45(@(z,Q) -Q^2*sqrt(2*(Q-log(Q)-1)), tspan, y0);
end

% Overlay the leading-order perturbation solution.
plot(-z,Q,'r')
shg

xlim([-5,0])
ylim([0,2])
xlabel('z')
ylabel('Density')
legend('PDE solution','Leading order perturbation')
title('PDE solution v leading order perturbation solution')
box on
shg


print(gcf,'-depsc2',[filepath_save_figs '\' 'qmin_evolution_ctm_late_perturb' num2str(t_loop) '.eps'])
saveas(gcf,[filepath_save_figs '\' 'qmin_evolution_ctm_late_perturb' num2str(t_loop) '.fig'])
saveas(gcf,[filepath_save_figs '\' 'qmin_evolution_ctm_late_perturb' num2str(t_loop) '.jpg'])

%% Phase planes

function_PLOTS_phase_plane(q_hist,L_hist,t_hist, kappa, phi,filepath_save_figs,simulation_id)


end