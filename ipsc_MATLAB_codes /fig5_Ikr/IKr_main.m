%Main File to generate Figure 5: IKr Models
close all
clear
Vm_plot=-150:50;

%% Wu Lab Dataset Specifc model
load Wu_ikr
[tau_act_Wu, act_inf_Wu]=gating_calculation([x_KR(2:5),x_KR(10)], Vm_plot);
[tau_inact_Wu, inact_inf_Wu]=gating_calculation([x_KR(6:9), x_KR(11)], Vm_plot);

%% Ma et al. Dataset Specifc model
load Ma_ikr
[tau_act_ma, act_inf_ma]=gating_calculation([x_KR(2:5),x_KR(10)], Vm_plot);
[tau_inact_ma, inact_inf_ma]=gating_calculation([x_KR(6:9), x_KR(11)], Vm_plot);

%% Es Salah Lamoureux et al. Dataset Specifc model
load EsSalahLamoureux_ikr
[tau_act_Essalah, act_inf_Essalah]=gating_calculation([x_KR(2:5),x_KR(10)], Vm_plot);
[tau_inact_Essalah, inact_inf_Essalah]=gating_calculation([x_KR(6:9), x_KR(11)], Vm_plot);

%% Bellin et al. Dataset Specifc model:
load Bellin_ikr
[tau_act_Bellin, act_inf_Bellin]=gating_calculation([x_KR(2:5),x_KR(10)], Vm_plot);
[tau_inact_Bellin, inact_inf_Bellin]=gating_calculation([x_KR(6:9), x_KR(11)], Vm_plot);

%% Baseline IKr Model
load Baseline_ikr
[tau_act_Baseline, act_inf_Baseline]=gating_calculation([x_KR(2:5),x_KR(10)], Vm_plot);
[tau_inact_Baseline, inact_inf_Baseline]=gating_calculation([x_KR(6:9), x_KR(11)], Vm_plot);

%% Plot Fig 5A: Steady State
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, act_inf_ma,'color', [0 .45 .74]);
plot(Vm_plot, act_inf_Wu, 'color', [.85 .33 .1]);
plot(Vm_plot, act_inf_Essalah, 'color', [.49 .18 .56]);
plot(Vm_plot, act_inf_Bellin, 'color', [.47 .67 .19]);
plot(Vm_plot, act_inf_Baseline, 'color', [0 0 0]);
plot(Vm_plot, inact_inf_ma,'color', [0 .45 .74]);
plot(Vm_plot, inact_inf_Wu, 'color', [.85 .33 .1]);
plot(Vm_plot, inact_inf_Essalah, 'color', [.49 .18 .56]);
plot(Vm_plot, inact_inf_Bellin, 'color', [.47 .67 .19]);
plot(Vm_plot, inact_inf_Baseline, 'color', [0 0 0]);
xlabel('Voltage (mV)');
ylabel('Normalized I_{Kr}')
legend('Ma et al.','Wu Lab', 'Es Salah Lamoureux et al.', 'Bellin et al.', 'Baseline')
legend boxoff

%% Plot Fig 5C: Tau activation
figure,set(gcf,'color','w')
plot(Vm_plot, tau_act_ma,'color', [0 .45 .74]);
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, tau_act_Wu, 'color', [.85 .33 .1]);
plot(Vm_plot, tau_act_Essalah, 'color', [.49 .18 .56]);
plot(Vm_plot, tau_act_Bellin, 'color', [.47 .67 .19]);
plot(Vm_plot, tau_act_Baseline, 'color', [0 0 0]);
xlabel('Voltage (mV)');
ylabel('Tau_{act} (ms)');
legend('Ma et al.','Wu Lab', 'Es Salah Lamoureux et al.', 'Bellin et al.', 'Baseline')
legend boxoff


%% Plot Fig 5D: Tau inactivation
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, tau_inact_ma,'color', [0 .45 .74]);
plot(Vm_plot, tau_inact_Wu, 'color', [.85 .33 .1]);
plot(Vm_plot, tau_inact_Essalah, 'color', [.49 .18 .56]);
plot(Vm_plot, tau_inact_Bellin, 'color', [.47 .67 .19]);
plot(Vm_plot, tau_inact_Baseline, 'color', [0 0 0]);
xlabel('Voltage (mV)');
ylabel('Tau_{inact} (ms)');
legend('Ma et al.','Wu Lab', 'Es Salah Lamoureux et al.', 'Bellin et al.', 'Baseline')
legend boxoff


%% Function to calculate gating variable properties
function [ tau_x, x_inf] = gating_calculation(  var, V)
x1=var(1); x2=var(2); x5=var(3); x6=var(4);
x4=1/((1/x2)+(1/x6));x3=x5*x1;

alpha=x1.*exp(V./x2); beta=x3.*exp(V./x4);
 
tau_x=(1./(alpha+beta))+var(5);
x_inf=alpha./(alpha+beta);
end
