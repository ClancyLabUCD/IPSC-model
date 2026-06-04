%Generates Fig12D: If
close all
clear all
clc

load Param_If
xF_spont=all_parameters_spont(:,77:82);
all_params_XF=all_parameters(:,77:82);
variations=length(xF_spont(:,1));

cell2plot=8000;
if cell2plot>variations
    cell2plot=variations;
end

Vm_plot=-150:2:0;

%% calculate kinetic properties for viable spontaneously beating cells
[tau_act, act_inf]=gating_calculation(xF_spont(:,2:6), Vm_plot);

%% calculate kinetic properties for ALL randomly determined parameterizations
[tau_act_all,act_inf_all ]=gating_calculation(all_params_XF(:,2:6), Vm_plot);

%% Plot for fig 12C Tau act
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, tau_act_all(1:cell2plot,:), 'color', [.5 .5 .5],  'LineWidth',.5);
plot(Vm_plot, tau_act(1:cell2plot,:), 'color', [0 .45 .74],  'LineWidth',.5);
ylabel('Tau_{act} (ms)')
xlabel('Voltage (mV)')

%% Plot for fig 12D Steady State
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, act_inf_all(1:cell2plot,:), 'color', [.5 .5 .5],  'LineWidth',.5);
plot(Vm_plot, act_inf(1:cell2plot,:), 'color', [0 .45 .74],  'LineWidth',.5);
ylabel('Normalized I_{F}')
xlabel('Voltage (mV)')


%% Function to calculate gating variable properties
function [ tau_x, x_inf] = gating_calculation(  var, V)
x1=var(:,1); x2=var(:,2); x5=var(:,3); x6=var(:,4);
x4=1./((1./x2)+(1./x6));x3=x5.*x1;

alpha=x1.*exp(V./x2); beta=x3.*exp(V./x4);
 
tau_x=(1./(alpha+beta))+var(:,5);
x_inf=alpha./(alpha+beta);
end


