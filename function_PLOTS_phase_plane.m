function function_PLOTS_phase_plane(q_hist,L_hist,t_hist, kappa, phi,filepath_save_figs,simulation_id)


%% PDE solution
% for forward travelling wave use the final time.

%obtain q and L
q_final_time =q_hist(:,end);
L_final_time = L_hist(end);

%calculate p

%dz from ctm discretisation
dz=0.00001; %for improved accuracy
%dz=0.0001; %for speed, solution has good accuracy for parameters considered
nodesz=round(1/dz) +1;


change_in_z = L_final_time/nodesz;
diff_q_final_time = (q_final_time(2:end)-q_final_time(1:end-1))/change_in_z;
p_final_time = (1./q_final_time(2:end).^2).*diff_q_final_time;

figure
plot(q_final_time(2:end),p_final_time,'c','LineWidth',4)

% for negative wavespeeds define plot_times_perturb as a time when the travelling wave solution exists (after transient time to form and before boundary effects)


%% Calculate wavespeed for a forward travelling wave

c = (L_hist(end)-L_hist(end-1))/(t_hist(end)-t_hist(end-1));

% Y(1) is Q
% Y(2) is p

% System of ODES given by Eq. (6) 
f = @(t,Y) [Y(2)*Y(1)^2;
    -Y(1)*( c*Y(2)*Y(1) + (1-Y(1)))];

% System of ODES given by Eq. (6) with time in reverse
fminus = @(t,Y) [-Y(2)*Y(1)^2;
    Y(1)*( c*Y(2)*Y(1) + (1-Y(1)))];

y1 = linspace(0,1.5,10);


%% plot the boundary condition
%p=-c*Q

hold on
plot(y1, -c*y1, 'g')

shg

%% plot the other boundary condition

hold on
plot(y1, (1/phi)*(1 - kappa*y1), 'b')

shg


%% overlay the trajectory associated with the travelling wave solution

% determine the initial condition based on linear stability analysis
if kappa < 1
    
    y1_init= 1 + 0.0000001*1;
    y2_init = 0 + 0.0000001*(2/(c+sqrt(c^2 + 4)));
    
    
elseif kappa == 1
    
    y1_init = 1;
    y2_init = 0;
    
elseif kappa > 1
    
    y1_init= 1 - 0.0000001*1;
    y2_init = 0 - 0.0000001*(2/(c+sqrt(c^2 + 4)));
    
end

% plot the trajectory
hold on
for y20 = y2_init
    [ts,ys] = ode45(f,[0,100],[y1_init;y2_init]);
    plot(ys(:,1),ys(:,2),'m')
    plot(ys(1,1),ys(1,2),'ko') % starting point
    plot(ys(end,1),ys(end,2),'ks') % ending point
end


%% overlay the trajectories corresponding to eigenvectors 

% using linear stability analysis determine the initial condition and then plot the corresponding trajectory.
if kappa < 1
    
      hold on
    y1_init= 1 - 0.0001*1;
    y2_init = 0 - 0.0001*(2/(c+sqrt(c^2 + 4)));
    
    for y20 = y2_init
        [ts,ys] = ode45(f,[0,6.5],[y1_init;y2_init]);
        plot(ys(:,1),ys(:,2),'k')
    end
    
    hold on
    y1_init= 1 - 0.0001*1;
    y2_init = 0 - 0.0001*(2/(c-sqrt(c^2 + 4)));
    
    for y20 = y2_init
        [ts,ys] = ode45(fminus,[0,8.4],[y1_init;y2_init]);
        plot(ys(:,1),ys(:,2),'k')
    end
    
    
    hold on
    y1_init= 1 + 0.0001*1;
    y2_init = 0 + 0.0001*(2/(c-sqrt(c^2 + 4)));
    
    for y20 = y2_init
        [ts,ys] = ode45(fminus,[0,8],[y1_init;y2_init]);
        plot(ys(:,1),ys(:,2),'k')
    end
    
   
    
elseif kappa > 1
    
     hold on
    y1_init= 1 + 0.0001*1;
    y2_init = 0 + 0.0001*(2/(c+sqrt(c^2 + 4)));
    
    for y20 = y2_init
        [ts,ys] = ode45(f,[0,9.5],[y1_init;y2_init]);
        plot(ys(:,1),ys(:,2),'k')
    end
    
    hold on
    y1_init= 1 - 0.0001*1;
    y2_init = 0 - 0.0001*(2/(c-sqrt(c^2 + 4)));
    
    for y20 = y2_init
        [ts,ys] = ode45(fminus,[0,6.2],[y1_init;y2_init]);
        plot(ys(:,1),ys(:,2),'k')
    end
    
    
    hold on
    y1_init= 1 + 0.0001*1;
    y2_init = 0 + 0.0001*(2/(c-sqrt(c^2 + 4)));
    
    for y20 = y2_init
        [ts,ys] = ode45(fminus,[0,5.7],[y1_init;y2_init]);
        plot(ys(:,1),ys(:,2),'k')
    end
    shg
    
    
end


yticks([-0.5,0,0.5])
ylim([-0.5,0.5])
xlim([0,1.5])
xlabel('Q')
ylabel('p')
title('(Q,p) phase plane')
legend('PDE Solution','Free boundary condition','Free boundary condition','Trajectory')
shg


%% save the figure

print(gcf,'-depsc2',[filepath_save_figs '\' 'Qp_phase_plane' simulation_id '.eps'])
saveas(gcf,[filepath_save_figs '\' 'Qp_phase_plane' simulation_id '.fig'])
saveas(gcf,[filepath_save_figs '\' 'Qp_phase_plane' simulation_id '.jpg'])



end
