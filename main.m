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

num_eng_stg1 = 5; % Number of engines each stage
num_eng_stg2 = 1;
num_eng_stg3 = 1;

% --- Hydrogen Engine
HY_Chamber_Pres = 150 * 1e5; % Pa, chamber pressure
HY_throat_area = 0.1; % m2, throat area
R_c_HY = 0.2; % m, engine chamber radius
L_c_HY = 0.5; % m, engine chamber axial length
trunc_length_HY = 15; % m, the length we truncate off the engine nozzle

% --- Hydrocarbon(LCH4) Engine
HC_Chamber_Pres = 150 * 1e5; % Pa
HC_throat_area = 0.1; % m2, throat area
R_c_HC = 0.2; % m, engine chamber radius
L_c_HC = 0.5; % m, engine chamber axial length
trunc_length_HC = 0.2; % m, the length we truncate off the engine nozzle



%-----------------Engine Analysis -----------------% 
% Using data from CEA, get things like R_specific, Chamber Temperature, Gamma (We can designate a pressure, say 350 bars) 
% Or we can adjust the chamber pressure until engine mass is reasonable
% --- Engine type 1 (Hydrogen, upper stage)
in.gamma = 1.1443; % data from CEA
in.T0 = 3731.9; % K, data from CEA
in.R_specific = 1000*R/14.480; % J K-1 kg-1, specific gas constant of exhaust

in.P0 = HY_Chamber_Pres; % Pa, hydrogen engine chamber pressure
in.At = HY_throat_area; % Throat area of hydrogen engine
in.throat_width = sqrt(HY_throat_area); % Assume square nozzle
[Me_HY, Te_HY, Pe_HY, ue_HY, Ae_At_HY, F_HY, Isp_HY, xw, yw, trunc_index, Prat, Xmesh] = calcEngine(in, 'vacuum', trunc_length_HY, scale);
in.Me = Me_HY;
m_engine_HY = calcEngineMass(in, trunc_index, rho_wall, R_c_HY, L_c_HY, sigma_h, xw, yw, scale, Prat, Xmesh);

% disp([Me_HY, Te_HY, Pe_HY, ue_HY, Ae_At_HY, F_HY, Isp_HY])

% --- Engine type 2 (Hydrocarbon - LCH4, lower 3 stages)
in.gamma = 1.1333; % data from CEA
in.T0 = 3739.3; % K, data from CEA
in.R_specific = 1000*R/24.560; % J K-1 kg-1, specific gas constant of exhaust

in.P0 = HC_Chamber_Pres; % Pa, hydrogen engine chamber pressure
in.At = HC_throat_area; % Throat area of hydrogen engine
in.throat_width = sqrt(HC_throat_area); % Assume square nozzle
[Me_HC, Te_HC, Pe_HC, ue_HC, Ae_At_HC, F_HC, Isp_HC, xw, yw, trunc_index, Prat, Xmesh]=calcEngine(in, 'sea level', trunc_length_HC, scale);
in.Me = Me_HC;
m_engine_HC = calcEngineMass(in, trunc_index, rho_wall, R_c_HC, L_c_HC, sigma_h, xw, yw, scale, Prat, Xmesh);
% disp([Me_HC, Te_HC, Pe_HC, ue_HC, Ae_At_HC, F_HC, Isp_HC])

%---------------Upper Stage Analysis---------------%
% Pavan already done the hand calc (or code?) on Overleaf doc, just need to copy paste the algorithm over here
% Upper Stage Design for Mission from LEO to Moon and Back
% clear; clc;
% Given Data
Isp = 450; % specific impulse of the upper stage in seconds (LOX/LH2)
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
        if (Isp == 450)  % hydrogen
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

% To calculate the mass of fuel for the upper stage
mf_up = m01-mb1;
%disp('The mass of fuel for upper stage in Kg')
%disp(mf_up)


%---------------Lower Stages Analysis--------------%
% Compute optimal delta-V & mass distribution across 3 stages
% Pavan have you done any of this part?
% Below are the inputs for the function launch_vehicle_optimal_solution
x0 = [1400, 2600, 4000, -1e-3]; % Intial guess for Velocities of three stages and alpha value
Isp1 = Isp_HC; % Stage 1 Isp. Changes based on the propellant used
Isp2 = Isp_HC; % Stage 1 Isp. Changes based on the propellant used
Isp3 = Isp_HC; % Stage 1 Isp. Changes based on the propellant used
m_upper = 95820; % The wet mass of the upper stage.

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




% --------------------Plotting---------------------%




