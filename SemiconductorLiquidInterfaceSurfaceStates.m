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
% Numerically simulates the band diagrams and capacitive properties of surface states at the CuGaS2 solid-   %
% liquid interface: assuming n-type GuGaS2, a high supporting electrolyte concentration, and an akaline      %
% liquid region.  Hence, changes to system parameters for alternate materials and electrolytes, beyond those %
% values provided, may require direct changes to the source code to obtain physically sound solutions.       %
%                                                                                                            %
% References                                                                                                 %
%   1. A. Iqbal, Md. S. Hossain and K. H. Bevan, "The role of relative rate constants in determining surface %
%      state phenomena at semiconductor–liquid interfaces", Phys. Chem. Chem. Phys. 18, 29466 (2016).        %
%   2. B. Miao, A. Iqbal, K. Sangare, B. Marsan and K. H. Bevan, "Interpreting Interfacial Semiconductor-    %
%      Liquid Capacitive Characteristics Impacted by Surface States: A Theoretical and Experimental Study    %
%      of CuGaS2" (2020). 										     %
%                                                                                                            %
% ========================================================================================================== %

close all;
clear all;

% ===========================================================================================================
% ==================....Parameters, General Variables, and Operator Definitions Start....====================

% -----------------------------------------------------------------------------------------------------------
% Input the fundamental parameters
freq = '501Hz' % String specifying the relevant frequency
par_fundamental=load(cat(2,'fundamental_parameters',freq,'.txt'));% Input the fundamental parameters
par_ss=load(cat(2,'Fit',freq,'.txt'));                       % Input the surface state distribution

k = 8.617e-5;                   % Boltzmann constant [eV/K]
eps0 = 8.85e-12;                % Permittivity of free space [F/m]
q = 1.602e-19;                  % Charge on an electron [coul]
m = 9.11e-31;                   % Electron mass (kg)
planck_constant = 6.63e-34;     % [J.s] to calculate the "effective" density of the conduction band and valence band, Nc and Nv
N_avg = 6.023e23;               % Avogadro's number (# of particle/ mole)

T = par_fundamental(1);                       % Temperature in Kelvin

Ef = par_fundamental(2);                      % Fermi level [eV]      
EG = par_fundamental(3);                      % Enegy band gap [eV]   
KS_solid = par_fundamental(4);                % Dielectric constant   
me = par_fundamental(5);                      % Electron effective mass/electron mass  
mh = par_fundamental(6);                      % Hole effective mass/electron mass      
NA = par_fundamental(7);                      % p-type: Acceptor density (#/m^3)
ND = par_fundamental(8);                      % n-type: Donor density (#/m^3)

KS_liquid = par_fundamental(9);               % Dielectric constant for water of the liquid  
Compensating_cation = par_fundamental(10);    % Supporting electrolyte [mol]  
Compensating_anion = par_fundamental(11);     % Supporting electrolyte [mol]  

% -----------------------------------------------------------------------------------------------------------
% Calculate semiconductor parameters
kT = k*T;                                % in eV. unit of energy
Nc = 2*(2*pi*me*m*kT*q/(planck_constant^2))^(3/2);   % Effective density of conduction band states [#/m^3]
Nv = 2*(2*pi*mh*m*kT*q/(planck_constant^2))^(3/2);   % Effective density of valence band states [#/m^3]
ni_solid = sqrt(Nc*Nv)*exp(-EG/kT/2);    % Intrinsic carrier density (#/m^3)
mn = me*m;                               % Electron effective mass
mp = mh*m;                               % Hole effective mass

% -----------------------------------------------------------------------------------------------------------
% Calculate Liquid parameters (Alkaline solution assumed)
% This model assumes as a strongly akaline solution.
OH_minus = Compensating_anion;     % The supporting anions are OH-. The alkaline solution.
pOH = -log10(OH_minus);
pH = 14-pOH;
H_plus = 10^(-pH);
E_ox_SHE = 0-1*2.3*kT*(pH);        % Level for oxidized species (H+) with respect to SHE. (Ref. Chp. 2, R. Krol & M. Gratzel. Photoelectrochemical Hydrogen Production, Springer, 2012)  
E_red_SHE = 1.23+1*2.3*kT*(-pH);   % Level for reduced species (OH-) with repect to SHE. (Ref. Chp. 2, R. Krol & M. Gratzel. Photoelectrochemical Hydrogen Production, Springer, 2012)  
SHE_to_vac = -4.44;                % SHE to absolute scale conversion constant. (Ref. Chp 2, R. Krol & M. Gratzel. Photoelectrochemical Hydrogen Production, Springer, 2012) 
E_ox = SHE_to_vac-E_ox_SHE;        % SHE measure converted to absolute scale   
E_red = SHE_to_vac-E_red_SHE;      % SHE measure converted to absolute scale      
pH_pzc = 7;                        % pH at the point of zero charge (Ref. Chp 2, R. Krol & M. Gratzel. Photoelectrochemical Hydrogen Production, Springer, 2012) 
V_helmholtz = 2.3*kT*(pH-pH_pzc);  % Helmboltz potential drop (Ref. Chp 2, R. Krol & M. Gratzel. Photoelectrochemical Hydrogen Production, Springer, 2012) 
E_ox = E_ox - V_helmholtz;         % Shift redox levels relative to the band edges by the fixed Helmholtz potential drop
E_red = E_red - V_helmholtz;       % Shift redox levels relative to the band edges by the fixed Helmholtz potential drop


% -----------------------------------------------------------------------------------------------------------
% Spatial grid defintion for the semiconductor and liquid, depends on system screening, all lengths in meters.
% This grid should be modified if the screening length of the liquid and/or semiconductor are altered.
% Note: the semiconductor and liquid grid spacing should be the same at the interface (x=0) to avoid numerical 
%       errors when utilizing a variable dielectric constant across the interface 
x_1a = linspace(-.12e-6,-.05e-6,51);    % Extended region towards the bulk  (part a)
x_1b = linspace(x_1a(end),-6e-9,401);   % Extended region towards the bulk  (part b)
x_2 = linspace(x_1b(end),0,801);        % Define grid point in the high electric field region of the semiconductor
x_solid = [x_1a x_1b(2:end) x_2(2:end)];% Grid points on the semiconductor side
x_3 = linspace(0,0.03e-7,401);          % Define grid point in the Helmholtz layer of the liquid
x_4 = linspace(x_3(end),0.18e-7,51);    % Define grid point in the bulk of the liquid 
x_liquid = [x_3 x_4(2:end)];            % Grid point at the liquid side
x = [x_solid x_liquid(2:end)];          % Grid point in total

Nliquid = length(x_liquid)-1;           % Number of grid points inside liquid
Nsolid = length(x_solid);               % Number of grid points inside solid

N = length(x);                          % Total number of grid point  
h = zeros(1,N+1);                         
hd = diff(x);                           % Grid spacing
h = [hd(1) hd hd(end)];                 % Utilized in the Laplacian operator

% -----------------------------------------------------------------------------------------------------------
% Variable dielectric constant distribution
x_point = 20;                          % Number of grids point to accommodate dielectric constant change 
KSmat = x*0;
KSmat(1:Nsolid) = KS_solid;            % Dielectric constant matrix of semiconductor 
KSmat(Nsolid+1:N) = KS_liquid;         % Dielectric constant matrix of liquid
deps = 2.5e-11;  % making sure that the dielectric constant changes smoothly at the interface to avoid numerical error. 2.5e-11 m results in a 2 Ang transtion region between dielectric constants 
feps = 1./(1+exp(-x/deps));            % Assumes that the solid dielectric constant is less than the liquid
KSmat = feps*(KS_liquid - KS_solid) + KS_solid;   
eps = KSmat*eps0;                       % Overall variable dielectric constant in F/m

% -----------------------------------------------------------------------------------------------------------
% Modified Laplacian for variable dielectric
La = (2./(h(2:N).*(h(2:N)+h(3:N+1)))); % Lower Laplacian diagonal                                         
Lb = (2./(h(1:N).*h(2:N+1)));          % Centre Laplacian diagonal
Lc = 2./(h(2:N).*(h(1:N-1)+h(2:N)));   % Upper Laplacian diagonal
laplacian_dirichlet1 = -diag(Lb) + diag(Lc,1) + diag(La,-1);
laplacian_dirichlet = diag(eps)*laplacian_dirichlet1;

% -----------------------------------------------------------------------------------------------------------
% Gradient Operator
Ga = [0 1./(h(3:N)) 0];
Gc = [0 1./(h(3:N))];
grad_op = -diag(Ga) + diag(Gc,1);       
grad_coeff = (grad_op*eps')';
Ga_eps = [0 grad_coeff(2:N-1).*(1./(h(3:N))) 0];
Gc_eps = [0 grad_coeff(2:N-1).*(1./(h(3:N)))];
grad_eps = -diag(Ga_eps) + diag(Gc_eps,1);

% -----------------------------------------------------------------------------------------------------------
% Semiconductor charge parameters
ni = zeros(1,Nsolid);         % The grid point in semiconductor excluding the surface state (surface state is at the semiconductor side)
ni(1:end) = ni_solid;         % The intrisic carrier density of semiconductor  
nuclear_static = x*0;         % Composed of Donors & Accpetors 
accep = x*0;
accep(1:Nsolid) = -NA;        % p-type electrode has higher accecptor concentration
donor = x*0;
donor(1:Nsolid) = 1*ND;       % n-type electrode has higher donor concentration
nuclear_static = accep + donor; % Doping concentration (all nuclear charges)
ss_charge = x*0;

% -----------------------------------------------------------------------------------------------------------
% System setup before contact (bc), band alignments and liquid redox energies wrt vacuum at 0 eV
% -> Future versions need a general first-order method for defining the liquid "Fermi level" at all pH values.
Ef_solid_bc = -Ef*ones(1,Nsolid);                         % Fermi level
Ec_solid_bc = Ef_solid_bc-kT*log(ND/Nc);                  % Conduction band
Ev_solid_bc = Ec_solid_bc-EG;                             % Valence band
Ei_solid_bc = (Ec_solid_bc+Ev_solid_bc)/2-kT*log(Nc/Nv)/2;% Intrisic Fermi level 
Vac_bc = zeros(1,N);                                      % Use vacuum as reference level
Eox_liquid_bc = E_ox*ones(1,Nliquid);                     % Level for reduced species (H+)  
Ered_liquid_bc = E_red*ones(1,Nliquid);                   % Level for oxidized species (OH-)
Eredox_liquid_bc = Ered_liquid_bc;                        % Usually it's a photo-anode
EA = abs(Ec_solid_bc(1));                                 % Electron affinity

% -----------------------------------------------------------------------------------------------------------
% Surface state placement. The surface state distribution is described by overlapping Gaussians.
Ess_offset=par_ss(1,:);    % The mean energy of each surface state with respect to the conduction band [eV]
sigma_Gausian=par_ss(2,:); % The standard deviation of the Gaussian distribution
Area=par_ss(3,:);          % The area of the surface state 
Lss=max(par_ss(4,:))      % The depth of the surface state, assume they all have the same value
Iss = floor(Lss/(x(Nsolid)-x(Nsolid-1))); % Index of the surface state distribution      
Lss_Eg=linspace(0,-EG,150); 
fss_density = zeros(length(Ess_offset),length(Lss_Eg));

for ii=1:length(Ess_offset)
  Dss_density(ii,:) = 1./sqrt(2.*pi.*sigma_Gausian(ii).^2).*exp(-(Lss_Eg-Ess_offset(ii)).^2./(2.*sigma_Gausian(ii)^2));
  fss_density0(ii,:) = (Dss_density(ii,1:end-1)+Dss_density(ii,2:end)).*abs(diff(Lss_Eg))./2; 
  fss_density(ii,1:end-1) = fss_density0(ii,:);
  fss_density(ii,end) = fss_density(ii,end-1);   % The probability density of surface states per energy
end

if(size(fss_density,1)== 1)
  ss_density_per_grid = (1./(Lss.*Area)).*fss_density;  % surface state density per energy (mono-surface state)
else
  for ii=1:length(Ess_offset)
    ss_density_per_grid(ii,:) = (1./(Lss.*Area(ii))).*fss_density(ii,:);
  end
end

ss_density_per_energy = sum(ss_density_per_grid);  % surface state density per energy (multi-surface state)
ss_density=zeros(length(Lss_Eg),length(x));

% ====================....Parameters, General Variables, and Operator Definitions End....====================
% ===========================================================================================================


% ===========================================================================================================
% ============================....Potential Profile & Bias Variables Start....===============================
dV = 1*kT;                            % Potential Bias Step
Vmax = 5*dV;                          % Maximum forward bias
Vmin = -60*dV;                        % Maxium range of the reverse bias
V = cat(2,[0:dV:Vmax],[0:-dV:Vmin]);  % Restart the loop at zero bias before iterating to forward or reverse bias
numV = length(V);                     % Number of voltage steps
beta = 1*pi/90;                       % Gummel loop mixing parameter

phi_all = zeros(numV,N);     % Potential profile at a given bias
Ei_all = zeros(numV,Nsolid); % Intrinsic level in band gap at a given bias
Fn_all = zeros(numV,Nsolid); % Electron quasi-Fermi level at a given bias
Fp_all = zeros(numV,Nsolid); % Hole quasi-Fermi level at a given bias
Ec_all = zeros(numV,Nsolid); % Conduction band profile at a given bias
Ev_all = zeros(numV,Nsolid); % Valence band profile at a given bias
n_all = zeros(numV,N);       % Electron density at a given bias
p_all = zeros(numV,N);       % Hole density at a given bias
ss_all = zeros(numV,N);       % Net surface state charge density at a given bias
electrode_charge = zeros(numV,N); % Electrode charge density at a given bias (includes ss charge)
liquid_charge = zeros(numV,N);    % Liquid charge density at a given bias
cation_all = zeros(numV,N);       % Cation density at a given bias
anion_all = zeros(numV,N);        % Anion densit at a given bias
electrode_total = zeros(1,numV);  % Electrode charge at a given bias (includes ss charge)
liquid_total = zeros(1,numV);     % Liquid charge at a given bias
ss_total = zeros(1,numV);         % Surface state (ss) charge at a given bias
dEf_bc = Ef_solid_bc(1)-Eredox_liquid_bc(1) + V; % Fermi level offset between liquid and semiconductor bulk, 
                                                 % positive bias lowers the semiconductor Fermi level down 
                                                 % relative to the liquid Fermi level (and vice versa)
 
                                                  
L_sc = sqrt( ((2*KS_solid*eps0)/q) ...
            *(1/abs(nuclear_static(1)))*dEf_bc(1) ) % approximate semiconductor screening length at zero bias,
                                                    % see https://en.wikipedia.org/wiki/Depletion_region
   
Ngummel = 1000;  % Number of Gummel loop iterations in bias calculations

% =============================....Potential Profile & Bias Variables End....================================
% ===========================================================================================================


% ===========================================================================================================
% ==================....Potential Profile (Zero and Applied Bias) Calculations Start....=====================


% -----------------------------------------------------------------------------------------------------------
% Begin the bias loop
for iv=1:numV   %1:numV
   Vbias = V(iv);                                           % Applied bias                    
   printstring = sprintf('iv: %i | bias: %f', iv,Vbias);
   disp(printstring);
	   
   dEf = dEf_bc(iv);                                         % Offset between solid and liquid Fermi levels before contact
   Fn = (-EA-abs(kT*log(ND/Nc))-Vbias)*ones(1,Nsolid);       % Quasi-fermi level of electrons 
   Fp = (-EA-abs(kT*log(ND/Nc))-Vbias)*ones(1,Nsolid);       % Quasi-fermi level of holes  
   phi_sc_bulk = kT*log(abs(nuclear_static(1))/ni(1))-Fn(1); % Assumes an n-type electrode, computes step potential that would be needed to align Fermi levels with Ei as a reference
   phi_liq_bulk = phi_sc_bulk-dEf;                           % Potential at the bulk of the liquid
   if(abs(Vbias) < 1e-12)                 % Zero bias condition
     phi = zeros(1,N);                    % Electrostatic potential profile 
     phi(1:Nsolid) = phi_sc_bulk;         % Initial potential guess, solid side
     phi(Nsolid+1:N) = phi_liq_bulk;      % Initial potential guess, liquid side
   else                                   % Biased system, use previous guess and scale it
     phi = phi_all(iv-1,:)*abs(dEf_bc(iv)/dEf_bc(iv-1));
     phi = (phi-phi(1))+phi_sc_bulk;
     phi(end-2:end) = phi_liq_bulk;
   end
   phi_start = phi;

   % ------------------------------------------------------------------------------------------------------
   % Begin the Gummel loop
   for ii=1:1000 %Ngummel% Gummel loop start

    % .....................................................................................................
    % Semiconductor setup
      Ei = -phi(1:Nsolid);  % Assumes an n-type electrode        
      n = [ni.*exp((Fn-Ei)/kT) zeros(1,Nliquid)]; % electron concentration initial guess
      p = [ni.*exp((Ei-Fp)/kT) zeros(1,Nliquid)]; % hole concentration initial guess
      Ec = Ei-kT*log(ni./Nc);                     % conduction band initial guess
      Ev = Ei+kT*log(ni./Nv);                     % valence band initial guess
      nuclear = nuclear_static;
      solid_charge=p-n+nuclear;

    % .....................................................................................................
    % Surface state charge
      Ess = Ec(Nsolid-floor(Iss/2))+Lss_Eg; % Position of surface state
      fss = 1./(1+exp((Ess-Fn(1))/kT));
      ss_density_per_energy_new = ss_density_per_energy.*(1-fss); % Donor surface states assumed

      for j=1:length(Lss_Eg)
          for kk=Nsolid-Iss+1:Nsolid
            ss_density(j,kk) = ss_density_per_energy_new(j);  
          end
      end
      
      ss_charge = sum(ss_density);% surface charge density 
      solid_charge = solid_charge+ss_charge;
 
    % .....................................................................................................
    % Solution concentration perturbations from the bulk (currently assumes a strongly alkaline solution) !!!
      density_compensating_cation = [zeros(1,Nsolid) ...     % Supportiong cation concentration
                  Compensating_cation*1000*N_avg*exp(-(phi(Nsolid+1:N)-phi_liq_bulk)/kT)];     
      density_H_plus = [zeros(1,Nsolid) ...                  % Hydrogen (H+) ion concentration
                  H_plus*1000*N_avg*exp(-(phi(Nsolid+1:N)-phi_liq_bulk)/kT)];       
      cation = density_H_plus+density_compensating_cation;   % Cation's concentration                                                                           
      density_compensating_anion = [zeros(1,Nsolid) ...      % Supportiong anion concentration
                  Compensating_anion*1000*N_avg*exp((phi(Nsolid+1:N)-phi_liq_bulk)/kT)];        
      density_OH_minus = [zeros(1,Nsolid) ...                % Hydroxide (OH-) concentration
                  OH_minus*1000*N_avg*exp((phi(Nsolid+1:N)-phi_liq_bulk)/kT)];                            
      anion = density_OH_minus;           % Anion's concentration (assumes an alkaline solution) !!!

    % .....................................................................................................
    % LU decomposition solution to the potential
      A = laplacian_dirichlet+1*grad_eps - diag(q*(cation+anion+p+n));                                                                  
      A(1,:) = [1 zeros(1,N-1)];   % Use the Dirichlet boundary condition at the bulk of the sc and the bulk of the liquid
      A(N,:) = [zeros(1,N-1) 1];   % Use the Dirichlet boundary condition at the bulk of the sc and the bulk of the liquid
      vector = (-1*laplacian_dirichlet*(phi'))-1*grad_eps*(phi') + (q*(anion-cation-solid_charge))';
      vector(1) = 0;               % Dirichlet boundary condition continued
      vector(end) = 0;             % Dirichlet boundary condition continued
      dphi = (LU_dec(diag(A,-1),diag(A),diag(A,1),vector))';% LU decomposition 
      phi = phi+beta*(dphi');                                           
      tol = max(abs(dphi));                                         
      printstring = sprintf('-> ii: %i | tol: %f', ii,tol);
      if(tol<1e-5)% Test convergence
        break;
      end  

   end
   % End the Gummel loop
   % --------------------------------------------------------------------------------------------------------

   % --------------------------------------------------------------------------------------------------------
   % Store calculation results
   phi_all(iv,:) = phi;               % Potential
   Ei_all(iv,:) = -phi(1:Nsolid);     % Intrisic fermi level (The potential is defined as the negative intrisic Fermi level.)
   Fn_all(iv,:) = Fn;                 % Electron quasi-Fermi level
   Fp_all(iv,:) = Fp;                 % Hole quasi-Fermi level 
   Ec_all(iv,:) = Ei-(kT)*log(ni./Nc);% Conduction band
   Ev_all(iv,:) = Ei+(kT)*log(ni./Nv);% Valence band

   p_all(iv,:)=[ni.*exp((-Fp+Ei)/kT) zeros(1,Nliquid)];         % Hole concentration
   n_all(iv,:)=[ni.*exp((-Ei+Fn)/kT) zeros(1,Nliquid)];         % Electron concentration
   ss_all(iv,:)=ss_charge;
   electrode_charge(iv,:) = p_all(iv,:)-n_all(iv,:)+nuclear+ss_all(iv,:);% Charge on the electrode
   cation_all(iv,:) = cation;         % Cation concentration
   anion_all(iv,:) = anion;           % Anion concentration
   liquid_charge(iv,:) = cation-anion;% Charge on the electrode
   
 
   electrode_total(iv)=q*sum((electrode_charge(iv,1:end-1)+electrode_charge(iv,2:end))*.5.*diff(x));
   ss_total(iv)=q*sum((ss_charge(1:end-1)+ss_charge(2:end))*.5.*diff(x));
   liquid_total(iv)=q*sum((liquid_charge(iv,1:end-1)+liquid_charge(iv,2:end))*.5.*diff(x));
  
end

% ===================....Potential Profile (Zero and Applied Bias) Calculations End....======================
% ===========================================================================================================



% -----------------------------------------------------------------------------------------------------------
% Save the system data & plot it
filename = cat(2,'SemiconductorLiquidInterfaceSurfaceStates_',freq,'.mat');
save(filename);
PlottingDataSurfaceStates(freq);
