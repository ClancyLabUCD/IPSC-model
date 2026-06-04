%Main File to generate Figure 4: ICaL Models
close all
clear
Vm_plot=-100.1:50.1;

%% Ma et al. Dataset Specifc model
load Ma_ical
[tau_act_ma, act_inf_ma]=gating_calculation([x_cal(2:5),  x_cal(10)], Vm_plot);
[tau_inact_ma, inact_inf_ma]=gating_calculation([x_cal(6:9), x_cal(11)], Vm_plot);

%% Veerman et al. Dataset #1 Dataset Specifc model
load Veerman1_ical
[tau_act_veerman1, act_inf_veerman1]=gating_calculation([x_cal(2:5), x_cal(10)], Vm_plot);
[tau_inact_veerman1, inact_inf_veerman1]=gating_calculation([x_cal(6:9), x_cal(11)], Vm_plot);

%% Veerman et al. Dataset #2 Dataset Specifc model
load Veerman2_ical
[tau_act_veerman2, act_inf_veerman2] =gating_calculation([x_cal(2:5), x_cal(10)], Vm_plot);
[tau_inact_veerman2, inact_inf_veerman2]=gating_calculation([x_cal(6:9), x_cal(11)], Vm_plot);

%% Es Salah Lamoureux et al. Steady State Data:
load EsSalahLamoureux_ical
[tau_act_essalah, act_inf_essalah]=gating_calculation([x_cal(2:5), x_cal(10)], Vm_plot);
[tau_inact_essalah, inact_inf_essalah]=gating_calculation([x_cal(6:9), x_cal(11)], Vm_plot);

%% Baseline ICaL Model
load Baseline_ical
[tau_act_Baseline, act_inf_Baseline]=gating_calculation([x_cal(2:5), x_cal(10)], Vm_plot);
[tau_inact_Baseline, inact_inf_Baseline]=gating_calculation([x_cal(6:9), x_cal(11)], Vm_plot);

%% Plot Fig 4A: Steady State
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, act_inf_ma,'color', [0 .45 .74]);
plot(Vm_plot, act_inf_veerman1, 'color', [.85 .33 .1]);
plot(Vm_plot, act_inf_veerman2, 'color', [.49 .18 .56]);
plot(Vm_plot, act_inf_essalah, 'color', [.47 .67 .19]);
plot(Vm_plot, act_inf_Baseline, 'color', [0 0 0]);
plot(Vm_plot, inact_inf_ma, 'color', [0 .45 .74]);
plot(Vm_plot, inact_inf_veerman1, 'color', [.85 .33 .1]);
plot(Vm_plot, inact_inf_veerman2, 'color', [.49 .18 .56]);
plot(Vm_plot, inact_inf_essalah, 'color', [.47 .67 .19]);
plot(Vm_plot, inact_inf_Baseline, 'color', [0 0 0]);
xlabel('Voltage (mV)');
ylabel('Normalized I_{CaL}')
legend('Ma et al.','Veerman et al. (1)', 'Veerman et al. (2)','Es-Salah-Lamoureaux et al.', 'Baseline')
legend boxoff

%% Plot Fig 4C: Tau activation
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot,tau_act_ma,'color', [0 .45 .74]);
plot(Vm_plot,tau_act_veerman1, 'color', [.85 .33 .1]);
plot(Vm_plot,tau_act_veerman2, 'color', [.49 .18 .56]);
plot(Vm_plot,tau_act_essalah, 'color', [.47 .67 .19]);
plot(Vm_plot,tau_act_Baseline, 'color', [0 0 0]);
xlabel('Voltage (mV)');
ylabel('Tau_{act} (ms)');
legend('Ma et al.','Veerman et al. (1)', 'Veerman et al. (2)','Es-Salah-Lamoureaux et al.', 'Baseline')
legend boxoff

%% Plot Fig 4D: Tau inactivation
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot,tau_inact_ma,'color', [0 .45 .74]);
plot(Vm_plot,tau_inact_veerman1, 'color', [.85 .33 .1]);
plot(Vm_plot,tau_inact_veerman2, 'color', [.49 .18 .56]);
plot(Vm_plot,tau_inact_essalah, 'color', [.47 .67 .19]);
plot(Vm_plot,tau_inact_Baseline, 'color', [0 0 0]);
xlabel('Voltage (mV)');
ylabel('Tau_{inact} (ms)');
legend('Ma et al.','Veerman et al. (1)', 'Veerman et al. (2)','Es-Salah-Lamoureaux et al.', 'Baseline')
legend boxoff


%% Function to calculate gating variable properties
function [ tau_x, x_inf] = gating_calculation(  var, V)
x1=var(1); x2=var(2); x5=var(3); x6=var(4);
x4=1/((1/x2)+(1/x6));x3=x5*x1;

alpha=x1.*exp(V./x2); beta=x3.*exp(V./x4);
 
tau_x=(1./(alpha+beta))+var(5);
x_inf=alpha./(alpha+beta);
end

