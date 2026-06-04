%Generates Fig12B: ICaL
close all
clear all
clc

load Param_ICaL
xCaL_spont=all_parameters_spont(:,51:61);
all_params_XCaL=all_parameters(:,51:61);
variations=length(xCaL_spont(:,1));

cell2plot=8000;
if cell2plot>variations
    cell2plot=variations;
end

Vm_plot=-100:50;

%% calculate kinetic properties for viable spontaneously beating cells
[tau_act, act_inf ]=gating_calculation([xCaL_spont(:,2:5),xCaL_spont(:,10)], Vm_plot);
[tau_inact, inact_inf]=gating_calculation([xCaL_spont(:,6:9),xCaL_spont(:,11)], Vm_plot);

%% calculate kinetic properties for ALL randomly determined parameterizations
[tau_act_all, act_inf_all]=gating_calculation([all_params_XCaL(:,2:5),all_params_XCaL(:,10)], Vm_plot);
[tau_inact_all, inact_inf_all]=gating_calculation([all_params_XCaL(:,6:9),all_params_XCaL(:,11)], Vm_plot);

%% Plot for fig 12B Tau act
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, tau_act_all(1:cell2plot,:), 'color', [.5 .5 .5],  'LineWidth',.5);
plot(Vm_plot, tau_act(1:cell2plot,:), 'color', [0 .45 .74],  'LineWidth',.5);
xlabel('Voltage (mV)')
ylabel('Tau_{act} (ms)')

%% Plot for fig 12B Tau inact
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, tau_inact_all(1:cell2plot,:), 'color', [.5 .5 .5],  'LineWidth',.5);
plot(Vm_plot, tau_inact(1:cell2plot,:), 'color', [.47 .67 .19],  'LineWidth',.5);
ylabel('Tau_{inact} (ms)')
xlabel('Voltage (mV)')

%% Plot for fig 12B Steady State
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, act_inf_all(1:cell2plot,:), 'color', [.5 .5 .5],  'LineWidth',.5);
plot(Vm_plot, inact_inf_all(1:cell2plot,:), 'color', [.5 .5 .5],  'LineWidth',.5);
plot(Vm_plot, act_inf(1:cell2plot,:), 'color', [0 .45 .74],  'LineWidth',.5);
plot(Vm_plot, inact_inf(1:cell2plot,:), 'color', [.47 .67 .19],  'LineWidth',.5);
ylabel('Normalized I_{CaL}')
xlabel('Voltage (mV)')

%% Function to calculate gating variable properties
function [ tau_x, x_inf] = gating_calculation(  var, V)
x1=var(:,1); x2=var(:,2); x5=var(:,3); x6=var(:,4);
x4=1./((1./x2)+(1./x6));x3=x5.*x1;

alpha=x1.*exp(V./x2); beta=x3.*exp(V./x4);
 
tau_x=(1./(alpha+beta))+var(:,5);
x_inf=alpha./(alpha+beta);
end


