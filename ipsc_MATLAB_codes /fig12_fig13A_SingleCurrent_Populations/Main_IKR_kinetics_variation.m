%Generates Fig12C: IKr
close all
clear all
clc

load Param_IKr
xKR_spont=all_parameters_spont(:,23:33);
all_params_XKR=all_parameters(:,23:33);
variations=max(size(xKR_spont));

cell2plot=8000;
if cell2plot>variations
    cell2plot=variations;
end

Vm_plot=-150:50;

%% calculate kinetic properties for viable spontaneously beating cells
[tau_act,  act_inf]=gating_calculation([xKR_spont(:,2:5),xKR_spont(:,10)], Vm_plot);
[tau_inact, inact_inf]=gating_calculation([xKR_spont(:,6:9),xKR_spont(:,11)], Vm_plot);

%% calculate kinetic properties for ALL randomly determined parameterizations
[tau_act_all, act_inf_all ]=gating_calculation([all_params_XKR(:,2:5),all_params_XKR(:,10)], Vm_plot);
[tau_inact_all,  inact_inf_all]=gating_calculation([all_params_XKR(:,6:9),all_params_XKR(:,11)], Vm_plot);

%% Plot for fig 12C Tau act
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, tau_act_all(1:cell2plot,:), 'color', [.5 .5 .5],  'LineWidth',.5);
plot(Vm_plot, tau_act(1:cell2plot,:), 'color', [0 .45 .74],  'LineWidth',.5);
ylabel('Tau_{act} (ms)')
xlabel('Voltage (mV)')

%% Plot for fig 12C Tau inact
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, tau_inact_all(1:cell2plot,:), 'color', [.5 .5 .5],  'LineWidth',.5);
plot(Vm_plot, tau_inact(1:cell2plot,:), 'color', [.47 .67 .19],  'LineWidth',.5);
ylabel('Tau_{inact} (ms)')
xlabel('Voltage (mV)')

%% Plot for fig 12C Steady State
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, act_inf_all(1:cell2plot,:), 'color', [.5 .5 .5],  'LineWidth',.5);
plot(Vm_plot, inact_inf_all(1:cell2plot,:), 'color', [.5 .5 .5],  'LineWidth',.5);
plot(Vm_plot, act_inf(1:cell2plot,:), 'color', [0 .45 .74],  'LineWidth',.5);
plot(Vm_plot, inact_inf(1:cell2plot,:), 'color', [.47 .67 .19],  'LineWidth',.5);
ylabel('Normalized I_{KR}')
xlabel('Voltage (mV)')

%% Function to calculate gating variable properties
function [ tau_x, x_inf] = gating_calculation(  var, V)
x1=var(:,1); x2=var(:,2); x5=var(:,3); x6=var(:,4);
x4=1./((1./x2)+(1./x6));x3=x5.*x1;

alpha=x1.*exp(V./x2); beta=x3.*exp(V./x4);
 
tau_x=(1./(alpha+beta))+var(:,5);
x_inf=alpha./(alpha+beta);
end



