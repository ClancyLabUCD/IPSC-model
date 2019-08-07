% Kernik-Clancy iPSC-CM model
%**********************************************
%Kernik DC, Morotti S, Wu H, Garg P, Duff HJ, Kurokawa J, Jalife J, Wu JC, Grandi E, Clancy CE.
%A computational model of induced pluripotent stem-cell derived cardiomyocytes
%incorporating experimental variability from multiple data sources"
%J Physiol. 2019 Jul 6. doi: 10.1113/JP277724
%**********************************************
%
% Converted to C-code by Mao-Tsuen Jeng
%
% Colleen Clancy Lab @ UC davis
%
% May-21-2019
% Main File to generate Baseline Model Figures 10-11
%close all; clear; clc

Ys = load('ys.txt');
currents = load('currents.txt');

Time = Ys(:,1);
Cai = Ys(:,4);
Vm = Ys(:,2);

INaCa = currents(:,9);
IpCa = currents(:,10);
Iup = currents(:,15);

% %% Figure 10A & 10C: Calcium Flux analysis and Calcium Transient Trace
% ca_analysis_cpp( Time, Iup, INaCa, IpCa, Cai )
% 
% %% Figure 11A: action potential trace for baseline model 
% %figure,set(gcf,'color','w')
% figure(3);
% hold on;
% 
 plot(Time, Vm, 'b:', 'LineWidth', 2 ); %'Color', [.8 0 .18]);
 set(gca,'box','off','tickdir','out')
 ylabel('Voltage (mV)');
 xlabel('Time (ms)')
 
% hold off;


