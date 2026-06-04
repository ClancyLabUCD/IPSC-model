%Main file to generate Fig 15: Immature vs. Mature Analysis in 5-current
%variation population

clear
clc
close all

%% Input parameters
load ICs_5current
load Param_5current

immature_color=[0, .45, .74];
mature_color=[.64, .08, .18];

immature_color_all=[.3, .8, 1];
mature_color_all=[1, .6, .78];

%% Fig 15A: define immature and mature subpoulations
all_parameters_immature=zeros(size(all_parameters_spont));
all_parameters_mature=zeros(size(all_parameters_spont));
immature_count=0;
mature_count=0;

figure,set(gcf,'color','w')
hold on
for j=1:length(all_outputs_spont)
    if  all_outputs_spont(j,1)<75 && all_outputs_spont(j,4)<85 %immature
        immature_count=immature_count+1;
        plot(-1.*all_outputs_spont(j,1), all_outputs_spont(j,4), 'MarkerSize',5,'Marker','.','Color',immature_color )
        all_parameters_immature(immature_count,:)=all_parameters_spont(j,:);
    elseif all_outputs_spont(j,1)>75 && all_outputs_spont(j,4)>85 %mature
        mature_count=mature_count+1;
        plot(-1.* all_outputs_spont(j,1), all_outputs_spont(j,4), 'MarkerSize',5,'Marker','.','Color',mature_color )
        all_parameters_mature(mature_count,:)=all_parameters_spont(j,:);
    else
        plot(-1.*all_outputs_spont(j,1), all_outputs_spont(j,4), 'MarkerSize',5,'Marker','.','Color',[0.5 0.5 0.5] )
    end
end

all_parameters_immature=all_parameters_immature(1:immature_count,:);
all_parameters_mature=all_parameters_mature(1:mature_count,:);
set(gca,'box','off','tickdir','out')
ylabel('dv/dt (V/s)')
xlabel('MDP (mV)')
hold off

%% Fig 15C: Steady state inactivation curves for mature and immature subpopulations
Vm_plot=-150:50;
xNA_mature=all_parameters_mature(:,63:76);
xNA_immature=all_parameters_immature(:,63:76);

xNA_mature_avg=mean(xNA_mature);
xNA_immature_avg=mean(xNA_immature);

%INa mature AVERAGE:
[~, h_inf_mature_avg]=gating_calculation([xNA_mature_avg(:,6:9),xNA_mature_avg(:,13)], Vm_plot);
h_inf_mature_avg=h_inf_mature_avg.^2;

%INa immature AVERAGE:
[~, h_inf_immature_avg]=gating_calculation([xNA_immature_avg(:,6:9),xNA_immature_avg(:,13)], Vm_plot);
h_inf_immature_avg=h_inf_immature_avg.^2;

%INa mature Population:
[~, h_inf_mature]=gating_calculation([xNA_mature(:,6:9),xNA_mature(:,13)], Vm_plot);
h_inf_mature=h_inf_mature.^2;

%INa immature population:
[~, h_inf_immature]=gating_calculation([xNA_immature(:,6:9),xNA_immature(:,13)], Vm_plot);
h_inf_immature=h_inf_immature.^2;

figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, h_inf_immature, 'color', immature_color_all,  'LineWidth',.5);
plot(Vm_plot, h_inf_mature, 'color', mature_color_all,  'LineWidth',.5);
plot(Vm_plot, h_inf_mature_avg, 'color', mature_color,  'LineWidth',2);
plot(Vm_plot, h_inf_immature_avg, 'color', immature_color,  'LineWidth',2);
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


