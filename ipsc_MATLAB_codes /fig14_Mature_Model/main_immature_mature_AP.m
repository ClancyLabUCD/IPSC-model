%% Main File for iPSC_function.m
close all; clear; clc

immature_color=[0, .45, .74];
mature_color=[.64, .08, .18];
options = odeset('MaxStep',1,'InitialStep',2e-2);
run_time=1.4e3;
%% Run immature/baseline iPSC_function
load ICs_baseline
load baseline_parameter_inputs
[Time_immature, values_immature] = ode15s(@ipsc_function,[0, run_time],Y_init, options, baseline_parameter_inputs);
dvdt_immature=(values_immature(1:end-1,1)-values_immature(2:end,1))./(Time_immature(1:end-1)-Time_immature(2:end));
[~,dvdt_immature_max]=max(dvdt_immature);

%% Run mature iPSC_function
load ICs_mature
%implement experimental changes in maximal conductances for Ik1 and INa:
mature_parameter_inputs=baseline_parameter_inputs;
mature_parameter_inputs(17)=mature_parameter_inputs(17)*(11.24/5.67); %ik1
mature_parameter_inputs(63)=mature_parameter_inputs(63)*(187/129); %ina
mature_parameter_inputs(83)=1; %with Istim
%Run mature model:
[Time_mature, values_mature] = ode15s(@ipsc_function,[0, run_time],Y_init, options, mature_parameter_inputs);
dvdt_mature=(values_mature(1:end-1,1)-values_mature(2:end,1))./(Time_mature(1:end-1)-Time_mature(2:end));
[~,dvdt_mature_max]=max(dvdt_mature);

%% Figure 14A: action potential trace for immature and mature models
figure,set(gcf,'color','w')
hold on
plot(Time_immature-Time_immature(dvdt_immature_max), values_immature(:,1),'Color', immature_color);
plot(Time_mature-Time_mature(dvdt_mature_max), values_mature(:,1),'Color', mature_color);
set(gca,'box','off','tickdir','out')
ylabel('Voltage (mV)');
xlabel('Time (ms)')
legend('Immature','Mature')
legend boxoff


