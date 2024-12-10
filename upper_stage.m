% Upper Stage Design for Mission from LEO to Moon and Back
clear; clc;
% Given Data
Isp = 450; % specific impulse of the upper stage in seconds (LOX/LH2) Need to change it
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

disp('Payload Ratio for upper Satge:')
disp(PL_ratio)

% Mass of Payload + Cargo
PL_mass = 5000; %kg

% To calculate Intial Mass of the upper stage
PL_fraction = PL_ratio/(PL_ratio+1);
m01 = PL_mass/PL_fraction;

disp('Intial Mass of upper Stage / Wet Mass of upper stage (kg):')
disp(m01)

% To calculate Burn Out Mass or Dry Mass
mb1 = m01/MR_total;

disp('Burn Out Mass of Upper Stage / Dry Mass (kg):')
disp(mb1)

% To calculate the mass of fuel for the upper stage
mf_up = m01-mb1;
disp('The mass of fuel for upper stage in Kg')
disp(mf_up)