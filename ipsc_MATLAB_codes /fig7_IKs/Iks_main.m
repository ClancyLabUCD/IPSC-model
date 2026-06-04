%Main File to generate Figure 7: Iks Models
close all
clear
Vm_plot=-100:100;

%% Ma et al.Dataset Specifc model
load Ma_iks
[tau_act_Ma, act_inf_Ma]=gating_calculation(x_IKS(2:6), Vm_plot);
act_inf_Ma=act_inf_Ma.^2;

%% Ma, Wei et al. Patient Dataset Specifc model
load MaWei_patient_iks
[tau_act_MaWei_patient, act_inf_MaWei_patient]=gating_calculation(x_IKS(2:6), Vm_plot);
act_inf_MaWei_patient=act_inf_MaWei_patient.^2;

%% Ma, Wei et al. iCell Dataset Specifc model
load MaWei_icell_iks
[tau_act_MaWei_icell, act_inf_MaWei_icell]=gating_calculation(x_IKS(2:6), Vm_plot);
act_inf_MaWei_icell=act_inf_MaWei_icell.^2;

%% Baseline model of Iks
load Baseline_iks
[tau_act_Baseline, act_inf_Baseline]=gating_calculation(x_IKS(2:6), Vm_plot);
act_inf_Baseline=act_inf_Baseline.^2;

%% Plot Fig 7A: Steady State
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, act_inf_Ma,'color', [0 .45 .74]);
plot(Vm_plot, act_inf_MaWei_patient, 'color', [.85 .33 .1]);
plot(Vm_plot, act_inf_MaWei_icell, 'color', [.49 .18 .56]);
plot(Vm_plot, act_inf_Baseline, 'color', [0 0 0]);
xlabel('Voltage (mV)');
ylabel('Normalized I_{Ks}')
legend('Ma et al.','Ma, Wei et al. Patient', 'Ma, Wei et al. iCell', 'Baseline')
legend boxoff

%% Plot Fig 7B: Tau activation
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, tau_act_Ma,'color', [0 .45 .74]);
plot(Vm_plot, tau_act_MaWei_patient, 'color', [.85 .33 .1]);
plot(Vm_plot, tau_act_MaWei_icell, 'color', [.49 .18 .56]);
plot(Vm_plot, tau_act_Baseline, 'color', [0 0 0]);
xlabel('Voltage (mV)');
ylabel('tau_{act} (ms)');
legend('Ma et al.','Ma, Wei et al. Patient', 'Ma, Wei et al. iCell', 'Baseline')
legend boxoff

%% Function to calculate gating variable properties
function [ tau_x, x_inf] = gating_calculation(  var, V)
x1=var(1); x2=var(2); x5=var(3); x6=var(4);
x4=1/((1/x2)+(1/x6));x3=x5*x1;

alpha=x1.*exp(V./x2); beta=x3.*exp(V./x4);
 
tau_x=(1./(alpha+beta))+var(5);
x_inf=alpha./(alpha+beta);
end



