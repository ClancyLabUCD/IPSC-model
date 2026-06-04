% Main File to generate Baseline Model Figures 10-11
close all; clear; clc

load ICs_baseline
load baseline_parameter_inputs

%% Run iPSC_function
options = odeset('MaxStep',1,'InitialStep',2e-2);
run_time=3e3;
[Time, values] = ode15s(@ipsc_function,[0, run_time],Y_init, options, baseline_parameter_inputs);
Cai=values(:,3);
Vm=values(:,1);

%% Calculate select current traces:
INaCa = zeros(size(Time));
IpCa = zeros(size(Time));
Iup = zeros(size(Time));
for i= 1:size(values,1)
    [~, update_step_i] =  ipsc_function(Time(i), values(i,:),  baseline_parameter_inputs);    
    INaCa(i) = update_step_i(8);
    IpCa(i) = update_step_i(9);
    Iup(i) = update_step_i(14);
end

%% Figure 10A & 10C: Calcium Flux analysis and Calcium Transient Trace
ca_analysis( Time, Iup, INaCa, IpCa, Cai )

%% Figure 11A: action potential trace for baseline model 
figure,set(gcf,'color','w')
plot(Time, Vm,'Color', [.8 0 .18]);
set(gca,'box','off','tickdir','out')
ylabel('Voltage (mV)');
xlabel('Time (ms)')


