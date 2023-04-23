% ========================================================================================================== %
% This software was written in 2018 through 2020                                                             %
% by Asif Iqbal <asif.iqbal@mail.mcgill.ca>, Botong Maio <botong.miao@mail.mcgill.ca>                        %
% and Kirk H. Bevan <kirk.bevan@mcgill.ca> ("the authors").                                                  %
%                                                                                                            %
% LICENCE                                                                                                    %
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License.                    %
% This software is distributed without any warranty.                                                         %
% You should have received a copy of the Creative Commons Attribution-NonCommercial-ShareAlike 4.0           %
% International Public License along with this software.                                                     %
% If not, see <https://creativecommons.org/licenses/by-nc-sa/4.0/>.                                          %
%                                                                                                            %
% DISCLAIMER                                                                                                 %
% The authors and publishers make no warranties about the software, and disclaim liability for all uses of   %
% the software, to the fullest extent permitted by applicable law.  The authors do not recommend use of this %
% software for any purpose.  It is made freely available, solely to clarify points made in their analysis of %
% CuGaS2 as described below. When using or citing the software, you should not imply endorsement by the      %
% authors.                                                                                                   %
%                                                                                                            %
% DESCRIPTION                                                                                                %
% A separate function for plotting data obtained from "SemiconductorLiquidInterfaceSurfaceStates.m"          %
% ========================================================================================================== %

function [] = PlottingDataSurfaceStates(freq);

%%clear all;
%%close all;

% ----------------------------------------------------------------------------------------------------------
% Load the relevant data
%%freq = '50Hz';   % String specifying the relevant frequency
load(cat(2,'SemiconductorLiquidInterfaceSurfaceStates_',freq,'.mat')); % Input the simulation data
Exp_data = load(cat(2,'capacitanceEXP',freq,'.txt'));  % Input the experimental data

% ----------------------------------------------------------------------------------------------------------
% Dielectric constant plotting
figure; 
plot(x*1e9,KSmat,'x-'); 
xlabel('Postion (nm)');
ylabel('Relative Dielectric Constant');
title(' Relative Dielectric Constant Across the Solid-Liquid Interface');
grid on;

% ----------------------------------------------------------------------------------------------------------
% Plotting the DOS of the surface states relative to vacuum at the surface
figure;
plot(ss_density_per_energy.*1e-6, Ec_solid_bc(1)+Lss_Eg,... % [cm^-3]
     linspace(0,max(ss_density_per_energy.*1e-6),100),Ec_solid_bc(1)*ones(1,100),'k-',...
     linspace(0,max(ss_density_per_energy.*1e-6),100),Ev_solid_bc(1)*ones(1,100),'k-'); hold on;
s = cat(2,'Active Density Distribution of Surface States at ',freq);
title(s)
legend('DOS','E_C','E_V','location','northeastoutside');
xlabel('Density Distribution (cm^{-3})')
ylabel('Energy Relative to Vacuum at the Surface (eV)');

% ----------------------------------------------------------------------------------------------------------
% Plotting the bandstructure before contact (Position is in unit nm).
Eref_Ag = -4.44-0.21;   % redox potential of Ag/AgCl relative to 
figure;
plot(x*1e9,Vac_bc,'b',...
     x(1:Nsolid)*1e9, Ec_solid_bc,'k',...
     x(1:Nsolid)*1e9, Ef_solid_bc,'--k',...
     x(1:Nsolid)*1e9, Ev_solid_bc,'k',...
     x(1:Nsolid)*1e9, Ei_solid_bc,'--m',...
     x(1+Nsolid:N)*1e9, (Eredox_liquid_bc+V_helmholtz),'--r',...
     x(1+Nsolid:N)*1e9, (Eref_Ag).*ones(N-(1+Nsolid)+1,1) ,'g');
s = cat(2,'Energy Band Diagram with ',freq, ' DOS: Before Contact');
title(s);
legend('Vacuum','E_C','E_F','E_V','E_i','E_{F,L}','E_{Ag/AgCl}','location','northeastoutside');
xlabel('Position (nm)'); 
ylabel('Energy (eV)');
ylim([1.1*Ev_solid_bc(1) 0.1]);
set(gcf,'Position',[100 500 1000 500])
grid on;

% ----------------------------------------------------------------------------------------------------------
% Plotting the bandstructure after contact (equilibrium) (Position is in units of nm.)
figure;
phi_bd = -phi_all(1,:);    % band diagram (bd) electrostatic potential shifts
phi_bd = phi_bd - phi_bd(1) + (EA + Ec_all(1,1));
dhelholtz = 2.5000e-11;
phi_helmholtz = -V_helmholtz*(1./(1+exp(-(x-5*dhelholtz)/dhelholtz)));  % Place the helmholtz drop/rise just after the semiconductor region
plot(x*1e9, phi_bd+phi_helmholtz,'b',...
     x(1:Nsolid)*1e9, Ec_all(1,:),'k',...
     x(1:Nsolid)*1e9, Fn_all(1,:),'--k',...
     x(1:Nsolid)*1e9, Ev_all(1,:),'k',...
     x(1:Nsolid)*1e9, Ei_all(1,:),'--m',...
     x(1+Nsolid:N)*1e9, Eredox_liquid_bc+phi_bd(end),'--r',...
     x(1+Nsolid:N)*1e9, Eref_Ag.*ones(N-(1+Nsolid)+1,1)+(phi_bd(end)+phi_helmholtz(end)) ,'g'); hold on;
V_wrt_refAg = (Eref_Ag+phi_bd(end)+phi_helmholtz(end))-Fn_all(1,1);
%s = sprintf('Equilibrium Energy Band Diagram: 0 V External Bias at %3.3f V vs. Ag/AgCl',V_wrt_refAg);
s = sprintf(' 0 V External Bias at %3.3f V vs. Ag/AgCl',V_wrt_refAg);
s = cat(2,'Equilibrium Energy Band Diagram with ', freq, ' DOS:', s); 
title(s);
legend('Vacuum','E_C','E_F','E_V','E_i','E_{F,L}','E_{Ag/AgCl}','location','northeastoutside');
xlabel('Position (nm)'); 
ylabel('Energy (eV)')
ylim([1.1*Ev_solid_bc(1) 1.1*max(phi_bd)]);
set(gcf,'Position',[100 500 1000 500])
grid on;


% ----------------------------------------------------------------------------------------------------------
% Plotting the bandstructure under bias (Position is in units of nm).  % ip=55 is the bias of the flatband
ip = 35; %35; %56; %50; %66; %55; 
figure;
phi_bd_neq = -phi_all(ip,:);    % band diagram (bd) electrostatic potential shifts
phi_bd_neq = phi_bd_neq - phi_bd_neq(1) + (EA + Ec_all(ip,1));
zero_shift = phi_bd_neq(1);    % Shift to place O eV reference at the bulk electrode
plot(x*1e9, phi_bd_neq-zero_shift+phi_helmholtz,'b',...
     x(1:Nsolid)*1e9, Ec_all(ip,:)-zero_shift,'k',...
     x(1:Nsolid)*1e9, Fn_all(ip,:)-zero_shift,'--k',...
     x(1:Nsolid)*1e9, Ev_all(ip,:)-zero_shift,'k',...
     x(1:Nsolid)*1e9, Ei_all(ip,:)-zero_shift,'--m',...
     x(1+Nsolid:N)*1e9, Eredox_liquid_bc+phi_bd_neq(end)-zero_shift,'--r',...
     x(1+Nsolid:N)*1e9, Eref_Ag.*ones(N-(1+Nsolid)+1,1)+(phi_bd_neq(end)+phi_helmholtz(end))-zero_shift,'g'); hold on;
V_wrt_refAg = (Eref_Ag+phi_bd_neq(end)+phi_helmholtz(end))-Fn_all(ip,1);
s = sprintf(' %3.3f V External Bias at %3.3f V vs. Ag/AgCl',V(ip),V_wrt_refAg);
s = cat(2,'Biased Energy Band Diagram with ', freq, ' DOS:', s);
title(s);
legend('Vacuum','E_C','E_F','E_V','E_i','E_{F,L}','E_{Ag/AgCl}','location','northeastoutside');
xlabel('Position (nm)'); 
ylabel('Energy (eV)')
ylim([1.1*Ev_solid_bc(1) 1.1*max(phi_bd_neq)]);
set(gcf,'Position',[100 500 1000 500])
grid on;



% ----------------------------------------------------------------------------------------------------------
% Mott-Schottky Plot
phi_bd = -phi_all(1,:);    % band diagram (bd) electrostatic potential shifts
phi_bd = phi_bd - phi_bd(1) + (EA + Ec_all(1,1));
Fn_SHE = Fn_all(:,end) - (SHE_to_vac+phi_bd(end))+V_helmholtz;   % Offset of the electron quasi-Fermi level from the SHE reference electrode
    
[val,Vmax_index] = min(abs(V-Vmax));
[val,Vmin_index] = min(abs(V-Vmin));
Vindicies = cat(2,[Vmin_index:-1:Vmax_index+1],[2:Vmax_index]);
Vplot = V(Vindicies);
Fn_SHE_plot = -Fn_SHE(Vindicies);
Fn_SHE_mid = (Fn_SHE_plot(1:end-1)+Fn_SHE_plot(2:end))/2; 
Qsemi = electrode_total(Vindicies)-ss_total(Vindicies);
Csemi = (Qsemi(2:end)-Qsemi(1:end-1))./(Vplot(2:end)-Vplot(1:end-1));  % Differential capacitance
Qss = ss_total(Vindicies);
Css = (Qss(2:end)-Qss(1:end-1))./(Vplot(2:end)-Vplot(1:end-1));  % Differential capacitance
Qtot = electrode_total(Vindicies);
Ctot = (Qtot(2:end)-Qtot(1:end-1))./(Vplot(2:end)-Vplot(1:end-1));  % Differential capacitance
Qliquid = liquid_total(Vindicies);
Cliquid = (Qliquid(2:end)-Qliquid(1:end-1))./(Vplot(2:end)-Vplot(1:end-1)); 


% ----------------------------------------------------------------------------------------------------------
% The capacitance and Mott-Schottky plot with potential respect to Ag/AgCl reference
figure(6);
indicies = [11:length(Fn_SHE_mid)-5];  % Plotting indicies, to avoid huge capacitance in the semiconductor accumulation bias range
Fn_Ag = Fn_SHE_mid-0.21*ones(length(Fn_SHE_mid),1);  % Semiconductor electrochemical potential relative to Ag/AgCl
alpha = 14.12;   % Surface area scaling parameter
plot(Fn_Ag(indicies), alpha*(Ctot(indicies)*(100^(-2)))); hold on; % C/cm^2  
legend('Simulation result'); hold on;
title('Capacitance Plot');
xlabel('Potential (V vs Ag/AgCl)'); 
ylabel('Capacitance (Fcm^{-2})')
figure(7);
plot(Fn_Ag(indicies),1./((alpha*Ctot(indicies)*(100^(-2))).^2)); hold on; % C/cm^2 
legend('Simulation result','location','NorthWest'); hold on;
title('Mott-Schottky Plot');
xlabel('Potential (V vs Ag/AgCl)'); 
ylabel('1/C^2 (F^{-2}cm^{4})')



% ----------------------------------------------------------------------------------------------------------
% Plot experimental capacitance data atop simulation data
C_potential=Exp_data(1,:);
Csc_exp=Exp_data(2,:);
figure(6) 
plot(C_potential, Csc_exp,'g^'); hold on;
figure(7);
MS_exp=(1./Csc_exp).^2;
plot(C_potential, MS_exp,'g^'); hold on; 









