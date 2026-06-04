%Main file to generate Fig 13B-F: APs for 5-current variation populations
clear
clc
close all

%total cells ploted:
cell2plot=50; 
%to plot all cells from all populations (as shown in Fig. 13B), set all flags =1, and cell2plot=20000

%% Input parameters
load ICs_5current
load Param_5current

current5_color=[ .34 .88 .88];
%% Fig 13B: Run ode to generate APs and Plot
options = odeset('MaxStep',1,'InitialStep',2e-5);
figure,set(gcf,'color','w')
hold on
for j=1:cell2plot
    if j<length(all_ICs_spont)
        y0n =all_ICs_spont(j,:);
        model_parameter_inputs_trial = all_parameters_spont(j,:);
        tspan=[0 (60e3/(all_outputs_spont(j,5)).*2)-all_outputs_spont(j,3)]; %all_outputs_spont(j,5) = beating rate in bpm, %all_outputs_spont(j,3) = APD90
        [t,y] = ode15s(@ipsc_function,tspan,y0n,options,model_parameter_inputs_trial);
        [time_dvdtmax, endAP]= start_end_AP(t,y);
        plot(t(1:10:endAP)-time_dvdtmax,y(1:10:endAP,1), 'Color',current5_color *(.5 + .5.*rand));
    end     
end
set(gca,'box','off','tickdir','out')
ylabel('Voltage (mV)')
xlabel('Time (ms)')
hold off

%savefig('APs_5current.fig')

%% Fig 13D: Amp vs APD
figure,set(gcf,'color','w')
plot(all_outputs_spont(:,3), all_outputs_spont(:,2), 'MarkerSize',5,'Marker','.','LineStyle','none', 'Color',current5_color )
set(gca,'box','off','tickdir','out')
ylabel('Amplitude (mV)')
xlabel('APD_{90} (ms)')

%% Fig 13E: MDP vs APD
figure,set(gcf,'color','w')
plot(all_outputs_spont(:,3), -1.*all_outputs_spont(:,1), 'MarkerSize',5,'Marker','.','LineStyle','none', 'Color',current5_color )
set(gca,'box','off','tickdir','out')
ylabel('Max Diastolic Potential (mV)')
xlabel('APD_{90} (ms)')

%% Fig 13F: BPM vs dv/dt
figure,set(gcf,'color','w')
plot(all_outputs_spont(:,5), all_outputs_spont(:,4), 'MarkerSize',5,'Marker','.','LineStyle','none', 'Color',current5_color )
set(gca,'box','off','tickdir','out')
ylabel('Max Upstroke Velocity (V/s)')
xlabel('Beating Rate (bpm)')

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

