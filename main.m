%---------------------Constants----------------------%
g0 = 9.81;  % Gravity constant, m/s^2
R = 8.3145;  % Universal gas constant in J/(mol*K)



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






%-----------------Engine Analysis -----------------% 
% Using data from CEA, get things like R_specific, Chamber Temperature, Gamma (We can designate a pressure, say 350 bars) 
% Or we can adjust the chamber pressure until engine mass is reasonable
% --- Engine type 1 (Hydrogen)
Isp_HY = 1; %Placeholder

% --- Engine type 2 (Hydrocarbon - LCH4)
Isp_HC = 1; %Placeholder



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

