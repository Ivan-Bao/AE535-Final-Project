function [delta_ui, alpha, R1, R2, R3, lambda1, lambda2, lambda3, m_01, m_02, m_03] = optimize_dV_alpha(Isp1, fuel_type1, Isp2, fuel_type2, Isp3, fuel_type3, dV_total, m_upper) % Function definition for the optimization
    fun = @optimal; % function defined at the end
    
    % Initial guesses for the variables
    % x0 = [2666, 2666, 2667, -1e-3]; % Initial guess
    x0 = [dV_total/3, dV_total/3, dV_total/3, -1e-3];
    
    % Solve the optimization problem
    x = fsolve(fun, x0); % solve the optimization problem
    
    % Compute structural coefficients given optimized results.
    eps1 = structural_coefficient(fuel_type1, x(1));
    eps2 = structural_coefficient(fuel_type2, x(2));
    eps3 = structural_coefficient(fuel_type3, x(3));
    
    % Other parameters
    deltaV = dV_total; % total delta-v required
    g = 9.806; % acceleration due to gravity, m/s^2
    % m_upper = 21445; % total upper stage mass from the question 2
    
    % Mass ratios calculation
    R1 = (x(4)*Isp1*g + 1)/(x(4)*Isp1*g*eps1);
    R2 = (x(4)*Isp2*g + 1)/(x(4)*Isp2*g*eps2);
    R3 = (x(4)*Isp3*g + 1)/(x(4)*Isp3*g*eps3);
    
    % Calculate lambda values for the fuel mass fractions
    lambda1 = (1 - R1*eps1) / (R1 - 1);
    lambda2 = (1 - R2*eps2) / (R2 - 1);
    lambda3 = (1 - R3*eps3) / (R3 - 1);
    
    % Calculate the stage masses
    m_03 = m_upper * (1 + lambda3) / lambda3;
    m_02 = m_03 * (1 + lambda2) / lambda2;
    m_01 = m_02 * (1 + lambda1) / lambda1;
    
    % delta u for each stage
    delta_ui = x(1:3);  % Extracts delta u1, delta u2, delta u3 from the solution vector
    
    % Scaling factor
    alpha = x(4);  % Extracts alpha from the solution vector
    
    
    
    function F = optimal(x)
        % x(1): delta u1
        % x(2): delta u2
        % x(3): delta u3
        % x(4): alpha
    
        % params
        deltaV = 8200; % Total delta-v required
        g = 9.806; % Gravity, m/s^2
    
        % Isp1 = 350; % Isp for hydrocarbon (LOX/CH4)
        % Isp2 = 350; % Isp for hydrocarbon (LOX/CH4)
        % Isp3 = 450; % Isp for hydrogen (LOX/LH2)
    
        % Calculate structural coefficients epsilon for each stage
        if strcmp(fuel_type1, 'hydrogen')  % hydrogen
            eps1 = 0.07060945 + 0.1610852*exp(-0.849934*(0.001*x(1)));
        else  % hydrocarbon
            eps1 = 0.0305466 + 0.06892734*exp(-0.8586846*(0.001*x(1)));
        end
    
        if strcmp(fuel_type2, 'hydrogen')
            eps2 = 0.07060945 + 0.1610852*exp(-0.849934*(0.001*x(2)));
        else
            eps2 = 0.0305466 + 0.06892734*exp(-0.8586846*(0.001*x(2)));
        end
    
        if strcmp(fuel_type3, 'hydrogen')
            eps3 = 0.07060945 + 0.1610852*exp(-0.849934*(0.001*x(3)));
        else
            eps3 = 0.0305466 + 0.06892734*exp(-0.8586846*(0.001*x(3)));
        end
    
        % Calculate delta ui for each stage
        % delta ui = Isp*g*ln(Ri)
        % where Ri is the mass ratio, derived from the formula given
        F(1) = x(1) - Isp1*g*log((x(4)*Isp1*g + 1)/(alpha*Isp1*g*eps1));
        F(2) = x(2) - Isp2*g*log((x(4)*Isp2*g + 1)/(x(4)*Isp1*g*eps2));
        F(3) = x(3) - Isp3*g*log((x(4)*Isp3*g + 1)/(x(4)*Isp2*g*eps3));
    
        % Sum of each delta ui is deltaV
        F(4) = deltaV - x(1) - x(2) - x(3);
    
    end
end
