%---------------------Constants----------------------%
g0 = 9.81;  % Gravity constant, m/s^2
R = 8.3145;  % Universal gas constant in J/(mol*K)



%---------------------Parameters--------------------%
delta_v_maneuver = 1000 * [8.0, 0.2, 2.46, 1.48, 0.68, 0.14, 0.68, 0.68, 0.14, 3.14] ;
payload_mass = 1000; % kg
fuel_type_stg1 = 'hydrogen';


% Define parameters
fuel_type = 'hydrogen';
delta_v = 1000; % Example delta-v value in m/s

% Call the function
epsilon = utility_functions.structural_coefficient(fuel_type, delta_v)