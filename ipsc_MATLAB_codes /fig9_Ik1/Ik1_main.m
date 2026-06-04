% Main File to generate Figure 9: IK1 Models
close all
clear
Vm_plot=-150:50;

%% Kurokawa Lab Dataset Specifc model
load Kurokawa_ik1
[ IV_Ik1_kurokawa ] = ik1_IV(x_K1, Vm_plot,0);

%% Ma et al.Dataset Specifc model
load Ma_ik1
[ IV_Ik1_ma ] = ik1_IV(x_K1, Vm_plot,1);

%% Jalife Mature Dataset Specifc model
load Jalife_Mature_ik1
[ IV_Ik1_jalifemature ] = ik1_IV(x_K1, Vm_plot,2);

%% Jalife Immature Dataset Specifc model
load Jalife_Immature_ik1
[ IV_Ik1_jalifeimmature ] = ik1_IV(x_K1, Vm_plot,2);

%% Baseline Ik1 model
load Baseline_ik1
[ IV_Ik1_Baseline] = ik1_IV( x_K1, Vm_plot,2);

%% Plot IVs
figure,set(gcf,'color','w')
set(gca,'box','off','tickdir','out')
hold on
plot(Vm_plot, IV_Ik1_jalifemature, 'color', [.47 .67 .19]);
plot(Vm_plot, IV_Ik1_jalifeimmature, 'color', [.49 .18 .56]);
plot(Vm_plot,IV_Ik1_ma,'color', [0 .45 .74]);
plot(Vm_plot, IV_Ik1_kurokawa, 'color', [.85 .33 .1]);   
plot(Vm_plot, IV_Ik1_Baseline, 'color', [0 0 0]);
xlabel('Voltage (mV)');
ylabel('I_{K1} (pA/pF)');
legend(  'Jalife Mature', 'Jalife Immature','Ma et al.','Kurokawa Lab','Baseline Model')
legend boxoff

%% Function to calculate IV curve
function [ ik1 ] = ik1_IV(  var, Vm, protocol )

if protocol==1 %Ma
    Ki=150; %mM
elseif protocol==2 %Jalife
    Ki=148; %mM
else %Kurokawa
    Ki=125; %mM
end

Ko = 5.4;   % millimolar (in model_parameters)
R = 8.314472;   % joule_per_mole_kelvin (in model_parameters)
T = 310.0;   % kelvin (in model_parameters)'
F = 96.4853415;   % coulomb_per_mmole (in model_parameters)
E_K = R*T/F*log(Ko/Ki);

alpha=var(2).*exp((Vm+var(4))./var(3));beta=exp((Vm+var(6))./var(5));
x_inf=alpha./(alpha+beta);

g_k1=var(1);

ik1= g_k1.*sqrt(Ko/5.4).*x_inf.*(Vm-E_K);

end
