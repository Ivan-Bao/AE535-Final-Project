clear; clc; close all;

%---------------------Constants----------------------%
g0 = 9.81;  % Gravity constant, m/s^2
R = 8.3145;  % Universal gas constant in J/(mol*K)
sigma_h = 50 * 1e6; % Pa Chamber material stress limit, 50 MPa
rho_wall = 6 * 100^3 * 1e-3; % kg/m3 Chamber material density
scale = 10;


%---------------------Parameters--------------------%
% --- Upper stage
delta_v_up_stg = 1000 * [0.2, 2.46, 1.48, 0.68, 0.14, 0.68, 0.68, 0.14, 3.14];
payload_mass = 5000; % kg
delta_v_btm_three_stage = 1000 * 8.0;
num_eng_up_stg = 1; % For now assume upper stage has only one engine

% --- Bottom 3 stages
fuel_type_stg1 = 'hydrocarbon';
fuel_type_stg2 = 'hydrocarbon';
fuel_type_stg3 = 'hydrocarbon';
fuel_type_upper = 'hydrogen';

num_eng_stg1 = 5; % Number of engines each stage
num_eng_stg2 = 1;
num_eng_stg3 = 1;

% --- Hydrogen Engine
HY_Chamber_Pres = 160 * 1e5; % Pa, chamber pressure
HY_throat_area = 0.1; % m2, throat area
R_c_HY = 0.2; % m, engine chamber radius
L_c_HY = 0.5; % m, engine chamber axial length
trunc_length_HY = 20; % m, the length we truncate off the engine nozzle

% --- Hydrocarbon(LCH4) Engine
HC_Chamber_Pres = 160 * 1e5; % Pa
HC_throat_area = 0.2; % m2, throat area
R_c_HC = 0.2; % m, engine chamber radius
L_c_HC = 0.5; % m, engine chamber axial length
trunc_length_HC = 3; % m, the length we truncate off the engine nozzle

% --- RP1 Engine
RP_Chamber_Pres = 160 * 1e5; % Pa
RP_throat_area = 0.2; % m2, throat area
R_c_RP = 0.2; % m, engine chamber radius
L_c_RP = 0.5; % m, engine chamber axial length
trunc_length_RP = 3; % m, the length we truncate off the engine nozzle

% --- Hydrogen Sea level Engine
HYSea_Chamber_Pres = 160 * 1e5; % Pa, chamber pressure
HYSea_throat_area = 0.2; % m2, throat area
R_c_HYSea = 0.2; % m, engine chamber radius
L_c_HYSea = 0.5; % m, engine chamber axial length
trunc_length_HYSea = 3; % m, the length we truncate off the engine nozzle


%-----------------Engine Analysis -----------------% 
% Or we can adjust the chamber pressure until engine mass is reasonable
% --- Engine type 1 (Hydrogen, upper stage)

in.gamma = 1.1459; % data from CEA <----------------------------------------------------------------------------change these three with data at 150 bar chamber pressure
in.T0 = 3574.9; % K, data from CEA
in.R_specific = 1000*R/13.582; % J K-1 kg-1, specific gas constant of exhaust


in.P0 = HY_Chamber_Pres; % Pa, hydrogen engine chamber pressure
in.At = HY_throat_area; % Throat area of hydrogen engine
in.throat_width = sqrt(HY_throat_area); % Assume square nozzle
[Me_HY, Te_HY, Pe_HY, ue_HY, Ae_At_HY, mdot_HY, F_HY, Isp_HY, xw, yw, trunc_index, Prat, Xmesh, noz_len_1m, noz_ex_rad_1m] = calcEngine(in, 'vacuum', trunc_length_HY, scale);
in.Me = Me_HY;
m_engine_HY = calcEngineMass(in, trunc_index, rho_wall, R_c_HY, L_c_HY, sigma_h, xw, yw, scale, Prat, Xmesh);
nozzle_len_HY = sqrt(HY_throat_area/pi) * noz_len_1m; % Nozzle length
nozzle_dia_ex_HY = sqrt(HY_throat_area/pi) * noz_ex_rad_1m * 2;

% --- Engine type 2 (Hydrocarbon - LCH4, lower 3 stages)

in.gamma = 1.1307; % data from CEA <----------------------------------------------------------------------------change these three with data at 150 bar chamber pressure
in.T0 = 3676.7; % K, data from CEA
in.R_specific = 1000*R/23.36; % J K-1 kg-1, specific gas constant of exhaust


in.P0 = HC_Chamber_Pres; % Pa, methane engine chamber pressure
in.At = HC_throat_area; % Throat area of methane engine
in.throat_width = sqrt(HC_throat_area); % Assume square nozzle
[Me_HC, Te_HC, Pe_HC, ue_HC, Ae_At_HC, mdot_HC, F_HC, Isp_HC, xw, yw, trunc_index, Prat, Xmesh, noz_len_1m, noz_ex_rad_1m]=calcEngine(in, 'sea level', trunc_length_HC, scale);
in.Me = Me_HC;
m_engine_HC = calcEngineMass(in, trunc_index, rho_wall, R_c_HC, L_c_HC, sigma_h, xw, yw, scale, Prat, Xmesh);
nozzle_len_HC = sqrt(HC_throat_area/pi) * noz_len_1m; % Nozzle length
nozzle_dia_ex_HC = sqrt(HC_throat_area/pi) * noz_ex_rad_1m * 2;

% --- Engine type 3 (RP1, not used, for reference)


in.gamma = 1.1354; % data from CEA <----------------------------------------------------------------------------change these three to fit RP1, chamber pressure also at 150 bar
in.T0 = 3834.9; % K, data from CEA
in.R_specific = 1000*R/25.035; % J K-1 kg-1, specific gas constant of exhaust


in.P0 = RP_Chamber_Pres; % Pa, RP1 engine chamber pressure
in.At = RP_throat_area; % Throat area of RP1 engine
in.throat_width = sqrt(RP_throat_area); % Assume square nozzle
[Me_RP, Te_RP, Pe_RP, ue_RP, Ae_At_RP, mdot_RP, F_RP, Isp_RP, xw, yw, trunc_index, Prat, Xmesh, noz_len_1m, noz_ex_rad_1m]=calcEngine(in, 'sea level', trunc_length_RP, scale);
in.Me = Me_RP;
m_engine_RP = calcEngineMass(in, trunc_index, rho_wall, R_c_RP, L_c_RP, sigma_h, xw, yw, scale, Prat, Xmesh);
nozzle_len_RP = sqrt(RP_throat_area/pi) * noz_len_1m; % Nozzle length
nozzle_dia_ex_RP = sqrt(RP_throat_area/pi) * noz_ex_rad_1m * 2;

% --- Engine type 4 (Hydrgen, sea level optimized)
in.gamma = 1.1459; % data from CEA
in.T0 = 3574.9; % K, data from CEA
in.R_specific = 1000*R/13.582; % J K-1 kg-1, specific gas constant of exhaust

in.P0 = HYSea_Chamber_Pres; % Pa, Sea level hydrogen engine chamber pressure
in.At = HYSea_throat_area; % Throat area of sea level hydrogen engine
in.throat_width = sqrt(HYSea_throat_area); % Assume square nozzle
[Me_HYSea, Te_HYSea, Pe_HYSea, ue_HYSea, Ae_At_HYSea, mdot_HYSea ,F_HYSea, Isp_HYSea, xw, yw, trunc_index, Prat, Xmesh, noz_len_1m, noz_ex_rad_1m]=calcEngine(in, 'sea level', trunc_length_HYSea, scale);
in.Me = Me_HYSea;
m_engine_HYSea = calcEngineMass(in, trunc_index, rho_wall, R_c_HYSea, L_c_HYSea, sigma_h, xw, yw, scale, Prat, Xmesh);
nozzle_len_HYSea = sqrt(HYSea_throat_area/pi) * noz_len_1m; % Nozzle length
nozzle_dia_ex_HYSea = sqrt(HYSea_throat_area/pi) * noz_ex_rad_1m * 2;

%---------------Upper Stage Analysis---------------%
% Pavan already done the hand calc (or code?) on Overleaf doc, just need to copy paste the algorithm over here
% Upper Stage Design for Mission from LEO to Moon and Back
% clear; clc;
% Given Data
Isp = Isp_HY; % specific impulse of the upper stage in seconds (LOX/LH2)
g0 = 9.81; % gravitational constant in m/s^2

% Delta-V values for each maneuver in km/s (converted to m/s)
delta_Vs = [0.2, 2.46, 1.48, 0.68, 0.14, 0.68, 0.68, 0.14, 3.14] * 1000; % m/s
delta_Vup = sum(delta_Vs);

% Rocket equation to calculate mass ratio for each maneuver
MR = exp(delta_Vs ./ (Isp * g0));

% Calculate total mass ratio by multiplying all the individual mass ratios
MR_total = prod(MR);

% Display individual mass ratios and total mass ratio
disp('Total Mass Ratio for upper stage:');
disp(MR_total);


% Calculate structural coefficients epsilon for each stage
if (Isp > 320)  % hydrogen
    eps = 0.07060945 + 0.1610852*exp(-0.849934*(0.001*delta_Vup));
else  % hydrocarbon
    eps = 0.0305466 + 0.06892734*exp(-0.8586846*(0.001*delta_Vup));
end

disp('The structural coefficeint of upper stage :')
disp(eps)

% Calculating Payload Ratio
PL_ratio = (1-eps.* MR_total)./(MR_total - 1);

%disp('Payload Ratio for upper Satge:')
%disp(PL_ratio)

% Mass of Payload + Cargo
PL_mass = 5000; %kg

% To calculate Intial Mass of the upper stage
PL_fraction = PL_ratio/(PL_ratio+1);
m01 = PL_mass/PL_fraction;
%disp('Intial Mass of upper Stage / Wet Mass of upper stage (kg):')
%disp(m01)

% To calculate Burn Out Mass or Dry Mass
mb1 = m01/MR_total;
%disp('Burn Out Mass of Upper Stage / Dry Mass (kg):')
%disp(mb1)

% To calculate the mass of propellant for the upper stage
mp_up = m01-mb1;
%disp('The mass of propellant for upper stage in Kg')
%disp(mf_up)

% To calculate the structural mass of the upper stage
ms_up = mb1-PL_mass;
%disp('The structural mass of the upper stage in Kg:')
%disp(ms_up)

%---------------Lower Stages Analysis--------------%
% Compute optimal delta-V & mass distribution across 3 stages
% Pavan have you done any of this part?
% Below are the inputs for the function launch_vehicle_optimal_solution
x0 = [2667, 2667, 2667, -1e-3]; % Intial guess for Velocities of three stages and alpha value
Isp1 = Isp_HC; % Stage 1 Isp. Changes based on the propellant used
Isp2 = Isp_HC; % Stage 1 Isp. Changes based on the propellant used
Isp3 = Isp_HC; % Stage 1 Isp. Changes based on the propellant used
m_upper = 170950; % The wet mass of the upper stage.

% The function launch_vehicle_optimal_solution will return optimal delatV values for each stage, optimal alpha, structural coefficients, payload ratios, mass ratios, 
% wet mass of each stage, structural mass of each stage and propellant mass of each stgae.
[x_solution,eps1,eps2,eps3,lambda1,lambda2, lambda3, R1, R2, R3, m_1, m_2, m_3, m_s1, m_s2, m_s3, m_p1, m_p2, m_p3, delta_ui, alpha] =...
launch_vehicle_optimal_solution(x0, Isp1, Isp2, Isp3, m_upper);
% Number of engines per stages are defined earlier, with the following variable name.
% num_eng_stg1
% num_eng_stg2
% num_eng_stg3

% We probably need to check the thrust of each stage is reasonable? Since the bottom 3 stage is lift off from earth surface, 
% I figured the 8km/s delta V includes some cosine loss during lift off, but we still need to aim for a reasonable thrust-to-weight ratio until we reach orbit.
% Perhaps we can use F9 or Star ship's TWR as reference.


%% ----------------- Propellant Volume Calculations ------------------
% Densities of propellants (kg/m^3) from Sutton
reference_gravity = 10^3;

density_CH4 = 0.424 * reference_gravity;   % Liquid methane (CH4) [kg/m^3]
density_H2 = 0.071 * reference_gravity;   % Liquid hydrogen (H2)
density_O2 = 1.14 * reference_gravity;   % Liquid oxygen (LOX)

% Oxidizer-to-Fuel ratios (O/F)
OF_ratio_CH4 = 4.0;  % Typical O/F ratio for CH4
OF_ratio_H2 = 6.0;   % Typical O/F ratio for H2
OF_ratio_RP1 = 3;

% Stage 1 Propellant Volume
if strcmp(fuel_type_stg1, 'hydrocarbon')
    m_fuel1 = m_p1 / (1 + OF_ratio_CH4);
    m_oxidizer1 = m_p1 - m_fuel1;
    V_fuel1 = m_fuel1 / density_CH4;     % Fuel volume
    V_oxidizer1 = m_oxidizer1 / density_O2; % Oxidizer volume
    V_propellant1 = V_fuel1 + V_oxidizer1; % Total propellant volume
else % Hydrogen
    m_fuel1 = m_p1 / (1 + OF_ratio_H2);
    m_oxidizer1 = m_p1 - m_fuel1;
    V_fuel1 = m_fuel1 / density_H2;      % Fuel volume
    V_oxidizer1 = m_oxidizer1 / density_O2; % Oxidizer volume
    V_propellant1 = V_fuel1 + V_oxidizer1; % Total propellant volume
end

% Stage 2 Propellant Volume
if strcmp(fuel_type_stg2, 'hydrocarbon')
    m_fuel2 = m_p2 / (1 + OF_ratio_CH4);
    m_oxidizer2 = m_p2 - m_fuel2;
    V_fuel2 = m_fuel2 / density_CH4;
    V_oxidizer2 = m_oxidizer2 / density_O2;
    V_propellant2 = V_fuel2 + V_oxidizer2;
else % Hydrogen
    m_fuel2 = m_p2 / (1 + OF_ratio_H2);
    m_oxidizer2 = m_p2 - m_fuel2;
    V_fuel2 = m_fuel2 / density_H2;
    V_oxidizer2 = m_oxidizer2 / density_O2;
    V_propellant2 = V_fuel2 + V_oxidizer2;
end

% Stage 3 Propellant Volume
if strcmp(fuel_type_stg3, 'hydrocarbon')
    m_fuel3 = m_p3 / (1 + OF_ratio_CH4);
    m_oxidizer3 = m_p3 - m_fuel3;
    V_fuel3 = m_fuel3 / density_CH4;
    V_oxidizer3 = m_oxidizer3 / density_O2;
    V_propellant3 = V_fuel3 + V_oxidizer3;
else % Hydrogen
    m_fuel3 = m_p3 / (1 + OF_ratio_H2);
    m_oxidizer3 = m_p3 - m_fuel3;
    V_fuel3 = m_fuel3 / density_H2;
    V_oxidizer3 = m_oxidizer3 / density_O2;
    V_propellant3 = V_fuel3 + V_oxidizer3;
end

% Upper Stage Propellant Volume
m_p_upper = 154.22e3;
if strcmp(fuel_type_upper, 'hydrocarbon')
    m_fuel_upper = m_p_upper / (1 + OF_ratio_CH4);
    m_oxidizer_upper = m_p_upper - m_fuel_upper;
    V_fuel_upper = m_fuel_upper / density_CH4;
    V_oxidizer_upper = m_oxidizer_upper / density_O2;
    V_propellant_upper = V_fuel_upper + V_oxidizer_upper;
else % Hydrogen
    m_fuel_upper = m_p_upper / (1 + OF_ratio_H2);
    m_oxidizer_upper = m_p_upper - m_fuel_upper;
    V_fuel_upper = m_fuel_upper / density_H2;
    V_oxidizer_upper = m_oxidizer_upper / density_O2;
    V_propellant_upper = V_fuel_upper + V_oxidizer_upper;
end

% Display results
disp(" ")
disp('Stage 1 Propellant Volumes:');
disp(['Fuel Volume (m^3): ', num2str(V_fuel1)]);
disp(['Oxidizer Volume (m^3): ', num2str(V_oxidizer1)]);
disp(['Total Volume (m^3): ', num2str(V_propellant1)]);
disp(" ")

disp('Stage 2 Propellant Volumes:');
disp(['Fuel Volume (m^3): ', num2str(V_fuel2)]);
disp(['Oxidizer Volume (m^3): ', num2str(V_oxidizer2)]);
disp(['Total Volume (m^3): ', num2str(V_propellant2)]);
disp(" ")

disp('Stage 3 Propellant Volumes:');
disp(['Fuel Volume (m^3): ', num2str(V_fuel3)]);
disp(['Oxidizer Volume (m^3): ', num2str(V_oxidizer3)]);
disp(['Total Volume (m^3): ', num2str(V_propellant3)]);
disp(" ")

disp('Upper Stage Propellant Volumes:');
disp(['Fuel Volume (m^3): ', num2str(V_fuel_upper)]);
disp(['Oxidizer Volume (m^3): ', num2str(V_oxidizer_upper)]);
disp(['Total Volume (m^3): ', num2str(V_propellant_upper)]);

%%  ----------------- Injector Analysis: -------------------------------
% ------------ CH4 ------
% Injector Design with Dynamic Cd Selection
fprintf('\n====================\n');
fprintf('\t CH4 Engines Injector Analysis:\n\n');

% Define chamber conditions and propellant densities
P_c = HC_Chamber_Pres; % Chamber pressure in Pa
DeltaP_ratio = 0.2; % Pressure drop ratio
DeltaP = DeltaP_ratio * P_c; % Calculate pressure drop across injectors
rho_fuel = density_CH4; % Density of CH4 in kg/m^3
rho_oxidizer = density_O2; % Density of LOX in kg/m^3

% Total mass flow and O/F ratio
total_m_dot = 1.77e3; % Total mass flow rate in kg/s
O_F_ratio = OF_ratio_CH4; % Oxidizer-to-Fuel ratio for CH4 and LOX

% Calculate mass flow rates for fuel and oxidizer
m_dot_fuel = total_m_dot / (1 + O_F_ratio); % Calculate fuel mass flow rate
m_dot_oxidizer = total_m_dot - m_dot_fuel; % Remaining flow rate is oxidizer

% Injector orifice design parameters
orifice_type = 'short-tube-rounded'; % Specify orifice type (rounded short-tube)
C_d = 0.90; % Discharge coefficient for the selected orifice type
d_oxidizer_orifice = 3.18e-3; % Diameter of oxidizer orifice in meters

% Calculate oxidizer orifice area
A_oxidizer_orifice = pi * (d_oxidizer_orifice / 2)^2; % Area of a single oxidizer orifice in m^2

% Calculate oxidizer injection velocity and flow rate per orifice
v_oxidizer = C_d * sqrt(2 * DeltaP / rho_oxidizer); % Injection velocity for oxidizer in m/s
m_dot_oxidizer_per_orifice = C_d * A_oxidizer_orifice * sqrt(2 * rho_oxidizer * DeltaP); % Oxidizer flow rate per orifice

% Determine the required number of oxidizer injectors
num_oxidizer_injectors = ceil(m_dot_oxidizer / m_dot_oxidizer_per_orifice); % Total oxidizer injectors needed

% Match the number of fuel injectors to oxidizer injectors
num_fuel_injectors = num_oxidizer_injectors;

% Calculate required fuel orifice diameter to match injector count
m_dot_fuel_per_orifice = m_dot_fuel / num_fuel_injectors; % Fuel flow rate per orifice
A_fuel_orifice = m_dot_fuel_per_orifice / (C_d * sqrt(2 * rho_fuel * DeltaP)); % Area of a single fuel orifice in m^2
d_fuel_orifice = sqrt(4 * A_fuel_orifice / pi); % Diameter of fuel orifice in meters

% Display injector design results
fprintf('Selected orifice type: %s\n', orifice_type);
fprintf('Fuel injector orifice diameter: %.2f mm\n', d_fuel_orifice * 1e3);
fprintf('Oxidizer injector orifice diameter: %.2f mm\n', d_oxidizer_orifice * 1e3);
fprintf('Injection velocity (fuel): %.2f m/s\n', C_d * sqrt(2 * DeltaP / rho_fuel));
fprintf('Injection velocity (oxidizer): %.2f m/s\n', v_oxidizer);
fprintf('Number of fuel injectors: %d\n', num_fuel_injectors);
fprintf('Number of oxidizer injectors: %d\n', num_oxidizer_injectors);

% Chamber Design
% Define chamber throat parameters
A_throat = HC_throat_area; % Throat area in m^2
R_throat = sqrt(A_throat / pi); % Calculate throat radius from area

% Define and calculate chamber-to-throat area ratio
Ac_At = 1.3; % Chamber-to-throat area ratio
A_chamber = Ac_At * A_throat; % Calculate chamber cross-sectional area

% Calculate characteristic chamber length (L*) and actual chamber length
L_star = 80; % Characteristic chamber length in cm for LOX/CH4
L_chamber = (L_star / 100) * A_chamber / A_throat; % Actual chamber length in meters

% Display chamber design results
fprintf('Updated Characteristic chamber length (L*): %.2f cm\n', L_star);
fprintf('Updated Chamber area: %.4f m^2\n', A_chamber);
fprintf('Updated Chamber length: %.2f m\n', L_chamber);
fprintf('Updated Chamber to throat area ratio (Ac/At): %.2f\n', Ac_At);

% Add spacing between CH4 and H2 outputs for better console formatting
fprintf('\n====================\n\n');

% ------------ H2 ------
% Injector Design with Matching Injector Count
fprintf('\t H2 Engines Injector Analysis:\n\n');

% Define chamber conditions and propellant densities
P_c = HY_Chamber_Pres; % Chamber pressure in Pa for H2 engine
DeltaP_ratio = 0.2; % Pressure drop ratio
DeltaP = DeltaP_ratio * P_c; % Calculate pressure drop across injectors
rho_fuel = density_H2; % Density of H2 in kg/m^3
rho_oxidizer = density_O2; % Density of LOX in kg/m^3

% Total mass flow and O/F ratio
total_m_dot = 1.77e3; % Total mass flow rate in kg/s
O_F_ratio = OF_ratio_H2; % Oxidizer-to-Fuel ratio for LOX and H2

% Calculate mass flow rates for fuel and oxidizer
m_dot_fuel = total_m_dot / (1 + O_F_ratio); % Calculate fuel mass flow rate
m_dot_oxidizer = total_m_dot - m_dot_fuel; % Remaining flow rate is oxidizer

% Injector orifice design parameters
C_d = 0.90; % Discharge coefficient for the selected orifice type
d_oxidizer_orifice = 3.18e-3; % Diameter of oxidizer orifice in meters

% Calculate oxidizer orifice area
A_oxidizer_orifice = pi * (d_oxidizer_orifice / 2)^2; % Area of a single oxidizer orifice in m^2

% Calculate oxidizer injection velocity and flow rate per orifice
v_oxidizer = C_d * sqrt(2 * DeltaP / rho_oxidizer); % Injection velocity for oxidizer in m/s
m_dot_oxidizer_per_orifice = C_d * A_oxidizer_orifice * sqrt(2 * rho_oxidizer * DeltaP); % Oxidizer flow rate per orifice

% Determine the required number of oxidizer injectors
num_oxidizer_injectors = ceil(m_dot_oxidizer / m_dot_oxidizer_per_orifice); % Total oxidizer injectors needed

% Match the number of fuel injectors to oxidizer injectors
num_fuel_injectors = num_oxidizer_injectors;

% Calculate required fuel orifice diameter to match injector count
m_dot_fuel_per_orifice = m_dot_fuel / num_fuel_injectors; % Fuel flow rate per orifice
A_fuel_orifice = m_dot_fuel_per_orifice / (C_d * sqrt(2 * rho_fuel * DeltaP)); % Area of a single fuel orifice in m^2
d_fuel_orifice = sqrt(4 * A_fuel_orifice / pi); % Diameter of fuel orifice in meters

% Display injector design results
fprintf('Selected orifice type for H2: short-tube-rounded\n');
fprintf('Fuel injector orifice diameter: %.2f mm\n', d_fuel_orifice * 1e3);
fprintf('Oxidizer injector orifice diameter: %.2f mm\n', d_oxidizer_orifice * 1e3);
fprintf('Injection velocity (fuel): %.2f m/s\n', C_d * sqrt(2 * DeltaP / rho_fuel));
fprintf('Injection velocity (oxidizer): %.2f m/s\n', v_oxidizer);
fprintf('Number of fuel injectors: %d\n', num_fuel_injectors);
fprintf('Number of oxidizer injectors: %d\n', num_oxidizer_injectors);

% Chamber Design
% Define chamber throat parameters
A_throat = HY_throat_area; % Throat area in m^2 for H2 engine
R_throat = sqrt(A_throat / pi); % Calculate throat radius from area

% Define and calculate chamber-to-throat area ratio
Ac_At = 1.3; % Chamber-to-throat area ratio
A_chamber = Ac_At * A_throat; % Calculate chamber cross-sectional area

% Calculate characteristic chamber length (L*) and actual chamber length
L_star = 56; % Characteristic chamber length in cm for LOX/H2
L_chamber = (L_star / 100) * A_chamber / A_throat; % Actual chamber length in meters

% Display chamber design results
fprintf('Updated Characteristic chamber length (L*): %.2f cm\n', L_star);
fprintf('Updated Chamber area: %.4f m^2\n', A_chamber);
fprintf('Updated Chamber length: %.2f m\n', L_chamber);
fprintf('Updated Chamber to throat area ratio (Ac/At): %.2f\n', Ac_At);


% --------------------Plotting---------------------%




