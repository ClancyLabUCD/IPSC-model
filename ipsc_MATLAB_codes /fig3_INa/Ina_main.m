%Main File to generate Figure 3: Sodium Current Models
close all
clear
Vm_plot=-150:50;

%% Ma et al. Dataset Specifc model
load Ma_ina
[tau_m_ma, m_inf_ma ]=gating_calculation([x_NA(2:5), x_NA(12)], Vm_plot);
act_inf_ma=m_inf_ma.^3; %m^3
[tau_h_ma, h_inf_ma]=gating_calculation([x_NA(6:9), x_NA(13)], Vm_plot);
inact_inf_ma=h_inf_ma.^2; %(h_inf = j_inf)
tau_j_ma=gating_calculation([x_NA(10:11), x_NA(8:9), x_NA(14)], Vm_plot);

%% Jalife Lab Immature Dataset specific model
load Jalife_Immature_ina
[tau_m_Jalife_Immature, m_inf_Jalife_Immature]=gating_calculation([x_NA(2:5), x_NA(12)], Vm_plot);
act_inf_Jalife_Immature=m_inf_Jalife_Immature.^3; %m^3
[tau_h_Jalife_Immature, h_inf_Jalife_Immature]=gating_calculation([x_NA(6:9), x_NA(13)], Vm_plot);
inact_inf_Jalife_Immature=h_inf_Jalife_Immature.^2; %(h_inf = j_inf)
tau_j_Jalife_Immature=gating_calculation([x_NA(10:11), x_NA(8:9), x_NA(14)], Vm_plot);

%% Jalife Lab Mature Dataset specific model
load Jalife_Mature_ina
[tau_m_Jalife_Mature, m_inf_Jalife_Mature]=gating_calculation([x_NA(2:5), x_NA(12)], Vm_plot);
act_inf_Jalife_Mature=m_inf_Jalife_Mature.^3; %m^3
[tau_h_Jalife_Mature, h_inf_Jalife_Mature]=gating_calculation([x_NA(6:9), x_NA(13)], Vm_plot);
inact_inf_Jalife_Mature=h_inf_Jalife_Mature.^2; %(h_inf = j_inf)
tau_j_Jalife_Mature=gating_calculation([x_NA(10:11), x_NA(8:9), x_NA(14)], Vm_plot);

%% Baseline INa Model
load Baseline_ina
[tau_m_Baseline,m_inf_Baseline]=gating_calculation([x_NA(2:5), x_NA(12)], Vm_plot);
act_inf_Baseline=m_inf_Baseline.^3; %m^3
[tau_h_Baseline, h_inf_Baseline]=gating_calculation([x_NA(6:9), x_NA(13)], Vm_plot);
inact_inf_Baseline=h_inf_Baseline.^2; %(h_inf = j_inf)
tau_j_Baseline=gating_calculation([x_NA(10:11), x_NA(8:9), x_NA(14)], Vm_plot);

%% Plot Fig 3A: Steady State
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, act_inf_ma,'color', [0 .45 .74]);
plot(Vm_plot, act_inf_Jalife_Immature, 'color', [.49 .18 .56]);
plot(Vm_plot, act_inf_Jalife_Mature, 'color', [.47 .67 .19]);
plot(Vm_plot, act_inf_Baseline, 'color', [0 0 0]);
plot(Vm_plot, inact_inf_ma, 'color', [0 .45 .74]);
plot(Vm_plot, inact_inf_Jalife_Immature, 'color', [.49 .18 .56]);
plot(Vm_plot, inact_inf_Jalife_Mature, 'color', [.47 .67 .19]);
plot(Vm_plot, inact_inf_Baseline, 'color', [0 0 0]);
xlabel('Voltage (mV)');
ylabel('Normalized I_{Na}')
legend('Ma et al.', 'Jalife Immature','Jalife Mature', 'Baseline')
legend boxoff
hold off

%% Plot Fig 3C: Tau m
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, tau_m_ma,'color', [0 .45 .74]);
plot(Vm_plot, tau_m_Jalife_Immature, 'color', [.49 .18 .56]);
plot(Vm_plot, tau_m_Jalife_Mature, 'color', [.47 .67 .19]);
plot(Vm_plot, tau_m_Baseline, 'color', [0 0 0]);

xlabel('Voltage (mV)');
ylabel('Tau_m (ms)');
legend('Ma et al.', 'Jalife Immature','Jalife Mature', 'Baseline')
legend boxoff
hold off

%% Plot Fig 3D: Tau h
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, tau_h_ma,'color', [0 .45 .74]);
plot(Vm_plot, tau_h_Jalife_Immature, 'color', [.49 .18 .56]);
plot(Vm_plot, tau_h_Jalife_Mature, 'color', [.47 .67 .19]);
plot(Vm_plot, tau_h_Baseline, 'color', [0 0 0]);

xlabel('Voltage (mV)');
ylabel('Tau_h (ms)');
legend('Ma et al.', 'Jalife Immature','Jalife Mature', 'Baseline')
legend boxoff
hold off

%% Plot Fig 3E: Tau j
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, tau_j_ma,'color', [0 .45 .74]);
plot(Vm_plot, tau_j_Jalife_Immature, 'color', [.49 .18 .56]);
plot(Vm_plot, tau_j_Jalife_Mature, 'color', [.47 .67 .19]);
plot(Vm_plot, tau_j_Baseline, 'color', [0 0 0]);
xlabel('Voltage (mV)');
ylabel('Tau_j (ms)');
legend('Ma et al.', 'Jalife Immature','Jalife Mature', 'Baseline')
legend boxoff
hold off

%% Function to calculate gating variable properties
function [ tau_x, x_inf] = gating_calculation(  var, V)
x1=var(1); x2=var(2); x5=var(3); x6=var(4);
x4=1/((1/x2)+(1/x6));x3=x5*x1;

alpha=x1.*exp(V./x2); beta=x3.*exp(V./x4);
 
tau_x=(1./(alpha+beta))+var(5);
x_inf=alpha./(alpha+beta);
end



