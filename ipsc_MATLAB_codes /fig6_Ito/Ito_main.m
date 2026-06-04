%Main File to generate Figure 6: Ito Models
close all
clear
Vm_plot=-100:100;

%% Ma et al.Dataset Specifc model
load Ma_ito
[tau_act_Ma, act_inf_Ma]=gating_calculation([xTO(2:5),xTO(10)], Vm_plot);
[tau_inact_Ma, inact_inf_Ma]=gating_calculation([xTO(6:9),xTO(11)], Vm_plot);

%% Veerman et al.Dataset Specifc model
load Veerman_ito
[tau_act_Veerman, act_inf_Veerman]=gating_calculation([xTO(2:5),xTO(10)], Vm_plot);
[tau_inact_Veerman, inact_inf_Veerman]=gating_calculation([xTO(6:9),xTO(11)], Vm_plot);

%% Cordeiro et al.Dataset Specifc model
load Cordeiro_ito
[tau_act_Cordeiro, act_inf_Cordeiro]=gating_calculation([xTO(2:5),xTO(10)], Vm_plot);
[tau_inact_Cordeiro, inact_inf_Cordeiro]=gating_calculation([xTO(6:9),xTO(11)], Vm_plot);

%% Baseline Ito Model
load Baseline_ito
[tau_act_Baseline, act_inf_Baseline]=gating_calculation([xTO(2:5),xTO(10)], Vm_plot);
[tau_inact_Baseline, inact_inf_Baseline]=gating_calculation([xTO(6:9),xTO(11)], Vm_plot);

%% Plot Fig 6A: Steady State
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, act_inf_Ma,'color', [0 .45 .74] );
plot(Vm_plot, act_inf_Veerman, 'color', [.85 .33 .1] );
plot(Vm_plot, act_inf_Cordeiro, 'color',  [.47 .67 .19] );
plot(Vm_plot, act_inf_Baseline, 'color', [0 0 0] );
plot(Vm_plot, inact_inf_Ma,'color', [0 .45 .74] );
plot(Vm_plot, inact_inf_Veerman, 'color', [.85 .33 .1] );
plot(Vm_plot, inact_inf_Cordeiro, 'color',  [.47 .67 .19] );
plot(Vm_plot, inact_inf_Baseline, 'color', [0 0 0] );
xlabel('Voltage (mV)');
ylabel('Normalized I_{to}')
legend('Ma et al.', 'Veerman et al.','Cordeiro et al.', 'Baseline')
legend boxoff

%% Plot Fig 6C: Tau activation
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, tau_act_Ma,'color', [0 .45 .74] );
plot(Vm_plot, tau_act_Veerman, 'color', [.85 .33 .1] );
plot(Vm_plot, tau_act_Cordeiro, 'color',  [.47 .67 .19] );
plot(Vm_plot, tau_act_Baseline, 'color', [0 0 0] );
xlabel('Voltage (mV)');
ylabel('Tau_{act} (ms)')
legend('Ma et al.', 'Veerman et al.','Cordeiro et al.', 'Baseline')
legend boxoff

%% Plot Fig 6D: Tau inactivation
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, tau_inact_Ma,'color', [0 .45 .74] );
plot(Vm_plot, tau_inact_Veerman, 'color', [.85 .33 .1] );
plot(Vm_plot, tau_inact_Cordeiro, 'color',  [.47 .67 .19] );
plot(Vm_plot, tau_inact_Baseline, 'color', [0 0 0] );
xlabel('Voltage (mV)');
ylabel('Tau_{inact} (ms)')
legend('Ma et al.', 'Veerman et al.','Cordeiro et al.', 'Baseline')
legend boxoff

%% Function to calculate gating variable properties
function [ tau_x, x_inf] = gating_calculation(  var, V)
x1=var(1); x2=var(2); x5=var(3); x6=var(4);
x4=1/((1/x2)+(1/x6));x3=x5*x1;

alpha=x1.*exp(V./x2); beta=x3.*exp(V./x4);
 
tau_x=(1./(alpha+beta))+var(5);
x_inf=alpha./(alpha+beta);
end
