%Generates Fig12E: IK1
close all
clear all
clc

load Param_Ik1
xK1_spont=all_parameters_spont(:,17:22);
all_params_XK1=all_parameters(:,17:22);
variations=length(xK1_spont(:,1));

cell2plot=8000;
if cell2plot>variations
    cell2plot=variations;
end

Vm_plot=-150:100;

%% calculate kinetic properties for viable spontaneously beating cells
 [ ik1_plot ] = ik1_IV( xK1_spont, Vm_plot );

%% calculate kinetic properties for ALL randomly determined parameterizations
[ ik1_all_plot ] = ik1_IV( all_params_XK1, Vm_plot );

%% Plot for fig 12E IV 
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, ik1_all_plot(1:cell2plot,:), 'color', [.5 .5 .5],  'LineWidth',.5);
plot(Vm_plot, ik1_plot(1:cell2plot,:), 'color', [.85 .33 .1],  'LineWidth',.5);
ylabel('I_{k1} (pA/pF)')
xlabel('Voltage (mV)')

%% Function to calculate IV curve
function [ ik1 ] = ik1_IV(  var, Vm)
Ki=148;% millimolar
Ko = 5.4;   % millimolar (in model_parameters)
R = 8.314472;   % joule_per_mole_kelvin (in model_parameters)
T = 310.0;   % kelvin (in model_parameters)'
F = 96.4853415;   % coulomb_per_mmole (in model_parameters)
E_K = R*T/F*log(Ko/Ki);

alpha=var(:,2).*exp((Vm+var(:,4))./var(:,3));
beta=exp((Vm+var(:,6))./var(:,5));
x_inf=alpha./(alpha+beta);

g_k1=var(:,1);

ik1= g_k1.*sqrt(Ko/5.4).*x_inf.*(Vm-E_K);

end





