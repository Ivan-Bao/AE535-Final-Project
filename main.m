%---------------------Constants----------------------%
g0 = 9.81;  % Gravity constant, m/s^2
R = 8.3145;  % Universal gas constant in J/(mol*K)
sigma_h = 50 * 1e6; % Pa Chamber material stress limit, 50 MPa
rho_wall = 6 * 100^3 * 1e-3; % kg/m3 Chamber material density



%---------------------CEA Data----------------------%






%---------------------Parameters--------------------%
% --- Upper stage
delta_v_up_stg = 1000 * [0.2, 2.46, 1.48, 0.68, 0.14, 0.68, 0.68, 0.14, 3.14];
payload_mass = 5000; % kg
delta_v_btm_three_stage = 1000 * 8.0;
num_eng_up_stg = 1; % For now assume upper stage has only one engine

% --- Bottom 3 stages
fuel_type_stg1 = 'hydrogen';
fuel_type_stg2 = 'hydrocarbon';
fuel_type_stg3 = 'hydrocarbon';

num_eng_stg1 = 5; % Number of engines each stage
num_eng_stg2 = 1;
num_eng_stg3 = 1;

% --- Hydrogen Engine
HY_Chamber_Pres = 350 * 1e6; % Pa, chamber pressure
HY_throat_area = 0.05; % m2, throat area
R_c_HY = 0.5; % m, engine chamber radius
L_c_HY = 1; % m, engine chamber axial length

% --- Hydrocarbon(LCH4) Engine
HC_Chamber_Pres = 350 * 1e6; % Pa
HC_throat_area = 0.05; % m2, throat area
R_c_HC = 0.5; % m, engine chamber radius
L_c_HC = 1; % m, engine chamber axial length




%-----------------Engine Analysis -----------------% 
% Using data from CEA, get things like R_specific, Chamber Temperature, Gamma (We can designate a pressure, say 350 bars) 
% Or we can adjust the chamber pressure until engine mass is reasonable
% --- Engine type 1 (Hydrogen)
in.gamma = 1.4; % needed data from CEA
in.T0 = 3000; % K, needed data from CEA
in.R_specific = 200; % J K-1 kg-1, needed data from CEA (Computed from combustion molecular weight and universal gas constant)

in.P0 = HY_Chamber_Pres; % Pa, hydrogen engine chamber pressure
in.At = HY_throat_area; % Throat area of hydrogen engine
[Me_HY, Te_HY, Pe_HY, ue_HY, Ae_At_HY, F_HY, Isp_HY]=calcEngine(in, 'sea level'); % Might want some IF statements here to make sure the sea level and vacuum ones match the stage
m_engine_HY = calcEngineMass(in, rho_wall, R_c_HY, L_c_HY, sigma_h, 10, 0.1, 1);

% --- Engine type 2 (Hydrocarbon - LCH4)
in.gamma = 1.4; % needed data from CEA
in.T0 = 3000; % K, needed data from CEA
in.R_specific = 200; % J K-1 kg-1, needed data from CEA

in.P0 = HC_Chamber_Pres; % Pa, hydrogen engine chamber pressure
in.At = HC_throat_area; % Throat area of hydrogen engine
[Me_HC, Te_HC, Pe_HC, ue_HC, Ae_At_HC, F_HC, Isp_HC]=calcEngine(in, 'sea level');
m_engine_HC = calcEngineMass(in, rho_wall, R_c_HC, L_c_HC, sigma_h, 10, 0.1, 1);

% At this point all engine data including engine mass should be available

% After all engine numbers are obtained, assign Isp to each stage based on engine selection
if strcmp(fuel_type_stg1, 'hydrogen')
    Isp1 = Isp_HY;
else
    Isp1 = Isp_HC;
end

if strcmp(fuel_type_stg2, 'hydrogen')
    Isp2 = Isp_HY;
else
    Isp2 = Isp_HC;
end

if strcmp(fuel_type_stg3, 'hydrogen')
    Isp3 = Isp_HY;
else
    Isp3 = Isp_HC;
end
    


%---------------Upper Stage Analysis---------------%
m_upper = 1; % Placeholder, by the end of upper stage analysis we should have the upper stage mass
% Pavan already done the hand calc (or code?) on Overleaf doc, just need to copy paste the algorithm over here



%---------------Lower Stages Analysis--------------%
% Compute optimal delta-V & mass distribution across 3 stages
% Pavan have you done any of this part?
[delta_ui, alpha, R1, R2, R3, lambda1, lambda2, lambda3, m_01, m_02, m_03] = optimize_dV_alpha(Isp1, fuel_type1, Isp2, fuel_type2, Isp3, fuel_type3, delta_v_btm_three_stage, m_upper);
% Number of engines per stages are defined earlier, with the following variable name.
% num_eng_stg1
% num_eng_stg2
% num_eng_stg3

% We probably need to check the thrust of each stage is reasonable? Since the bottom 3 stage is lift off from earth surface, 
% I figured the 8km/s delta V includes some cosine loss during lift off, but we still need to aim for a reasonable thrust-to-weight ratio until we reach orbit.
% Perhaps we can use F9 or Star ship's TWR as reference.




% --------------------Plotting---------------------%




