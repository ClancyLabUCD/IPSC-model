% Main File to generate Figure 8: If Models
close all
clear
Vm_plot=-150.1:0.1;

%% Ma et al.Dataset Specifc model
load Ma_if
[tau_act_fit_Ma, act_inf_fit_Ma] =gating_calculation(x_F(2:6), Vm_plot);

%% Kurokawa lab Dataset Specifc model
load Kurokawa_if
[tau_act_fit_Kurokawa, act_inf_fit_Kurokawa] =gating_calculation(x_F(2:6), Vm_plot);

%% Baseline If model
load Baseline_if
[tau_act_fit_Baseline, act_inf_fit_Baseline] =gating_calculation(x_F(2:6), Vm_plot);

%% Plot Fig 8A: Steady State
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, act_inf_fit_Kurokawa,'color', [0 .45 .74]);
plot(Vm_plot, act_inf_fit_Ma,'color', [.85 .33 .1]);
plot(Vm_plot, act_inf_fit_Baseline,'color', [0 0 0]);
xlabel('Voltage (mV)');
ylabel('Normalized I_{f}')
legend('Kurokawa Lab', 'Ma et al.', 'Baseline')
legend boxoff

%% Plot Fig 8B: Tau activation
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot,tau_act_fit_Kurokawa,'color', [0 .45 .74]);
plot(Vm_plot,tau_act_fit_Ma,'color', [.85 .33 .1]);
plot(Vm_plot,tau_act_fit_Baseline,'color', [0 0 0]);
xlabel('Voltage (mV)');
ylabel('Tau_{act} (ms)');
legend('Kurokawa Lab', 'Ma et al.', 'Baseline')
legend boxoff

%% Function to calculate gating variable properties
function [ tau_x, x_inf] = gating_calculation(  var, V)
x1=var(1); x2=var(2); x5=var(3); x6=var(4);
x4=1/((1/x2)+(1/x6));x3=x5*x1;

alpha=x1.*exp(V./x2); beta=x3.*exp(V./x4);
 
tau_x=(1./(alpha+beta))+var(5);
x_inf=alpha./(alpha+beta);
end


