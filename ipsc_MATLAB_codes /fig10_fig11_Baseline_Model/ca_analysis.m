function [  ] = ca_analysis( Time, Iup, INaCa, IpCa, Ca )
%calcium transient analysis
%output: plots for Figure 10A&C

% Constants (copied from ipsc_function)
V_tot=3960; %um^3 from hwang et al.
Vc_tenT=16404; VSR_tenT=1094; V_tot_tenT=Vc_tenT+VSR_tenT;
Vc=V_tot*(Vc_tenT/V_tot_tenT);
Cm = 60; %pF
F = 96.4853415;   % coulomb_per_mmole (in model_parameters)

%% Find first beat to analyze
inds_time_800=find(Time>800);inds_time_800=inds_time_800(1);
inds_time_1600=find(Time>1600); inds_time_1600=inds_time_1600(1);
[~,inds1]=min(Ca(1:inds_time_800));
[~,inds2]=min(Ca(inds_time_800:inds_time_1600)); inds2=inds_time_800+inds2;

%% Calculate Normalized Ca2+ flux 
%take integral
intJserca = cumtrapz(Time,Iup);
intIncx_ca = cumtrapz(Time,-INaCa*2*Cm/(2.0*Vc*F));
intIpca = cumtrapz(Time,IpCa*Cm/(2.0*Vc*F));

%integral for first beat
fluxJserca = intJserca(inds1:inds2)-intJserca(inds1);
fluxIncx_ca = intIncx_ca(inds1:inds2)-intIncx_ca(inds1);
fluxIpca = intIpca(inds1:inds2)-intIpca(inds1);

%Normalize flux
flux_total=fluxJserca+fluxIncx_ca+fluxIpca;
ref=max(flux_total);
fluxJserca_norm=fluxJserca./ref;
fluxIncx_ca_norm=fluxIncx_ca./ref;
fluxIpca_norm=fluxIpca./ref;
Time_flux=Time(inds1:inds2);

%% plot figure 10A 
figure, set(gcf,'color','w')
plot((Time(inds1:inds2)-Time(inds1))./1000, Ca(inds1:inds2).*1e6, 'Color', [.8 0 .18])
set(gca,'box','off','tickdir','out')
legend boxoff
ylabel('[Ca^{2+}] (nM)')
xlabel('Time (s)');

%% plot figure 10C 
figure,set(gcf,'color','w')
plot((Time_flux-Time(inds1)),fluxJserca_norm,(Time_flux-Time(inds1)),fluxIncx_ca_norm,(Time_flux-Time(inds1)),fluxIpca_norm);
set(gca,'box','off','tickdir','out')
legend('SERCA', 'NCX', 'non-NCX')
legend boxoff
ylabel('Ca flux normalized')
xlabel('Time (ms)')

end

