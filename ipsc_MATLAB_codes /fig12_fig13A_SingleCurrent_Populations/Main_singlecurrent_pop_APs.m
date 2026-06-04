%Main file to generate Fig 12A: APs for single channel variation populations

clear
clc
close all

%total cells ploted:
cell2plot=200; 

%choose populations to plot:
plot_ik1=1;
plot_ikr=1;
plot_ical=1;
plot_ina=1;
plot_if=1;
%to plot all cells from all populations (as shown in Fig. 13A), set all flags =1, and cell2plot=40000

%% Input parameters
load ICs_If
load Param_If
If_ICs_spont= all_ICs_spont;
If_outputs_spont= all_outputs_spont;
If_parameters_spont= all_parameters_spont;

load ICs_ICaL
load Param_ICaL
ICaL_ICs_spont= all_ICs_spont;
ICaL_outputs_spont= all_outputs_spont;
ICaL_parameters_spont= all_parameters_spont;

load ICs_INa
load Param_INa
INa_ICs_spont= all_ICs_spont;
INa_outputs_spont= all_outputs_spont;
INa_parameters_spont= all_parameters_spont;

load ICs_IKr
load Param_IKr
IKr_ICs_spont= all_ICs_spont;
IKr_outputs_spont= all_outputs_spont;
IKr_parameters_spont= all_parameters_spont;

load ICs_Ik1
load Param_Ik1
IK1_ICs_spont= all_ICs_spont;
IK1_outputs_spont= all_outputs_spont;
IK1_parameters_spont= all_parameters_spont;

%% Run ode to generate APs and Plot
CaL_avg_color=[.49, .18, .56];
Na_avg_color=[0 .5 0];
Ik1_avg_color=[.8, .2, 0];
Ikr_avg_color=[.74 .53 0] ;
If_avg_color=[1 .1 1];

options = odeset('MaxStep',1,'InitialStep',2e-5);
figure,set(gcf,'color','w')
hold on
for j=1:cell2plot/(plot_ik1+plot_ikr+plot_ical+plot_ina+plot_if)
    
    %plot APs from If variation population
    if j<length(If_ICs_spont) && plot_if==1
        y0n =If_ICs_spont(j,:);
        model_parameter_inputs_trial = If_parameters_spont(j,:);
        tspan=[0 (60e3/(If_outputs_spont(j,5)).*2)-If_outputs_spont(j,3)]; %all_outputs_spont(j,5) = beating rate in bpm, %all_outputs_spont(j,3) = APD90
        [t,y] = ode15s(@ipsc_function,tspan,y0n,options,model_parameter_inputs_trial);
        [time_dvdtmax, endAP]= start_end_AP(t,y);
        plot(t(1:10:endAP)-time_dvdtmax,y(1:10:endAP,1), 'Color', If_avg_color);
    end    
    
    %plot APs from ICaL variation population
    if j<length(ICaL_ICs_spont) && plot_ical==1
        y0n =ICaL_ICs_spont(j,:);
        model_parameter_inputs_trial = ICaL_parameters_spont(j,:);
        tspan=[0 (60e3/(ICaL_outputs_spont(j,5)).*2)-ICaL_outputs_spont(j,3)]; %all_outputs_spont(j,5) = beating rate in bpm, %all_outputs_spont(j,3) = APD90
        [t,y] = ode15s(@ipsc_function,tspan,y0n,options,model_parameter_inputs_trial);
        [time_dvdtmax, endAP]= start_end_AP(t,y);
        plot(t(1:10:endAP)-time_dvdtmax,y(1:10:endAP,1), 'Color', CaL_avg_color);
    end
    
    %plot APs from INa variation population
    if j<length(INa_ICs_spont) && plot_ina==1
        y0n =INa_ICs_spont(j,:);
        model_parameter_inputs_trial = INa_parameters_spont(j,:);
        tspan=[0 (60e3/(INa_outputs_spont(j,5)).*2)-INa_outputs_spont(j,3)]; %all_outputs_spont(j,5) = beating rate in bpm, %all_outputs_spont(j,3) = APD90
        [t,y] = ode15s(@ipsc_function,tspan,y0n,options,model_parameter_inputs_trial);
        [time_dvdtmax, endAP]= start_end_AP(t,y);
        plot(t(1:10:endAP)-time_dvdtmax,y(1:10:endAP,1), 'Color', Na_avg_color);
    end
    
    %plot APs from IKr variation population
    if j<length(IKr_ICs_spont)
        y0n =IKr_ICs_spont(j,:);
        model_parameter_inputs_trial = IKr_parameters_spont(j,:);
        tspan=[0 (60e3/(IKr_outputs_spont(j,5)).*2)-IKr_outputs_spont(j,3)]; %all_outputs_spont(j,5) = beating rate in bpm, %all_outputs_spont(j,3) = APD90
        [t,y] = ode15s(@ipsc_function,tspan,y0n,options,model_parameter_inputs_trial);
        [time_dvdtmax, endAP]= start_end_AP(t,y);
        plot(t(1:10:endAP)-time_dvdtmax,y(1:10:endAP,1), 'Color', Ikr_avg_color);
    end
    
    %plot APs from IK1 variation population
    if j<length(IK1_ICs_spont) && plot_ik1==1
        y0n =IK1_ICs_spont(j,:);
        model_parameter_inputs_trial = IK1_parameters_spont(j,:);
        tspan=[0 (60e3/(IK1_outputs_spont(j,5)).*2)-IK1_outputs_spont(j,3)]; %all_outputs_spont(j,5) = beating rate in bpm, %all_outputs_spont(j,3) = APD90
        [t,y] = ode15s(@ipsc_function,tspan,y0n,options,model_parameter_inputs_trial);
        [time_dvdtmax, endAP]= start_end_AP(t,y);
        plot(t(1:10:endAP)-time_dvdtmax,y(1:10:endAP,1), 'Color', Ik1_avg_color);
    end
end
set(gca,'box','off','tickdir','out')
ylabel('Voltage (mV)')
xlabel('Time (ms)')
hold off

%savefig('APs_singlecurrentVar.fig')
%% function to determine range of AP to plot

function  [time_dvdtmax, endAP]= start_end_AP(t,y)
dvdt=(y(1:end-1,1)-y(2:end,1))./(t(1:end-1)-t(2:end));
[~,startAP]=max(dvdt(1:floor(.66*end)));
time_dvdtmax=t(startAP);
endAP=find(y(floor(.75*end):end,1)>-30);
if length(endAP)>1
    endAP=endAP(1)+floor(.75*length(y));
else
    endAP=length(y(:,1));
end
end

