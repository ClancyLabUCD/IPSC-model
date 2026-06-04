%Generates Fig12A: INa
close all
clear all

load Param_INa
xNA_spont=all_parameters_spont(:,63:76);
all_params_XNA=all_parameters(:,63:76);

variations=length(xNA_spont(:,1));

cell2plot=8000;
if cell2plot>variations
    cell2plot=variations;
end

Vm_plot=-150:50;

%% calculate kinetic properties for viable spontaneously beating cells
[tau_m, m_inf]=gating_calculation([xNA_spont(:,2:5),xNA_spont(:,12)], Vm_plot);
m_inf=m_inf.^3;
[tau_h,  h_inf ]=gating_calculation([xNA_spont(:,6:9),xNA_spont(:,13)], Vm_plot);
h_inf=h_inf.^2;
[tau_j]=gating_calculation([xNA_spont(:,10:11),xNA_spont(:,8:9), xNA_spont(:,14)], Vm_plot);

%% calculate kinetic properties for ALL randomly determined parameterizations
[tau_m_all,  m_inf_all]=gating_calculation([all_params_XNA(:,2:5),all_params_XNA(:,12)], Vm_plot);
 m_inf_all= m_inf_all.^3;
[tau_h_all,  h_inf_all]=gating_calculation([all_params_XNA(:,6:9),all_params_XNA(:,13)], Vm_plot);
 h_inf_all= h_inf_all.^2;
[tau_j_all]=gating_calculation([all_params_XNA(:,10:11),all_params_XNA(:,8:9), all_params_XNA(:,14)], Vm_plot);

%% Plot for fig 12A Tau m
figure,set(gcf,'color','w')
plot(Vm_plot, tau_m_all(1:cell2plot,:), 'color', [.5 .5 .5],  'LineWidth',.5);
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, tau_m_all(1:cell2plot,:), 'color', [.5 .5 .5],  'LineWidth',.5);
plot(Vm_plot, tau_m(1:cell2plot,:), 'color', [0 .45 .74],  'LineWidth',.5);
ylabel('Tau_{m} (ms)')
xlabel('Voltage (mV)')

%% Plot for fig 11A Tau h
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, tau_h_all(1:cell2plot,:), 'color', [.5 .5 .5],  'LineWidth',.5);
plot(Vm_plot, tau_h(1:cell2plot,:), 'color', [.47 .67 .19],  'LineWidth',.5);
ylabel('Tau_{h} (ms)')
xlabel('Voltage (mV)')

%% Plot for fig 12A Tau J
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, tau_j_all(1:cell2plot,:), 'color', [.5 .5 .5],  'LineWidth',.5);
plot(Vm_plot, tau_j(1:cell2plot,:), 'color', [.49 .18 .56],  'LineWidth',.5);
ylabel('Tau_{j} (ms)')
xlabel('Voltage (mV)')

%% Plot for fig 12A Steady State
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, m_inf_all(1:cell2plot,:), 'color', [.5 .5 .5],  'LineWidth',.5);
plot(Vm_plot, h_inf_all(1:cell2plot,:), 'color', [.5 .5 .5],  'LineWidth',.5);
plot(Vm_plot, m_inf(1:cell2plot,:), 'color', [0 .45 .74],  'LineWidth',.5);
plot(Vm_plot, h_inf(1:cell2plot,:), 'color', [.47 .67 .19],  'LineWidth',.5);
ylabel('Normalized I_{Na}')
xlabel('Voltage (mV)')

%% Function to calculate gating variable properties
function [ tau_x, x_inf] = gating_calculation(  var, V)
x1=var(:,1); x2=var(:,2); x5=var(:,3); x6=var(:,4);
x4=1./((1./x2)+(1./x6));x3=x5.*x1;

alpha=x1.*exp(V./x2); beta=x3.*exp(V./x4);
 
tau_x=(1./(alpha+beta))+var(:,5);
x_inf=alpha./(alpha+beta);
end




