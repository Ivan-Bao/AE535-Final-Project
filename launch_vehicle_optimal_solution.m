function [x_solution,eps1,eps2,eps3,lambda1,lambda2, lambda3, R1, R2, R3, m_1, m_2, m_3, m_s1, m_s2, m_s3, m_p1, m_p2, m_p3, delta_ui, alpha] = launch_vehicle_optimal_solution(x0, Isp1, Isp2, Isp3, m_upper)
    % x0: Initial guess for [delta u1, delta u2, delta u3, alpha]
    % Isp1, Isp2, Isp3: Specific impulses for the stages
    % m_upper: Mass of the upper stage

    % Solve the optimization problem using fsolve
    options = optimoptions('fsolve');
    x_solution = fsolve(@(x) optimal(x, Isp1, Isp2, Isp3), x0, options);

    % Extract delta u values and alpha
    delta_ui = x_solution(1:3);
    alpha = x_solution(4);

    % Calculate epsilon values
    eps1 = calc_eps(Isp1, delta_ui(1));
    eps2 = calc_eps(Isp2, delta_ui(2));
    eps3 = calc_eps(Isp3, delta_ui(3));

    % Constants
    g = 9.806; % Gravitational constant

    % Mass ratios
    R1 = ((alpha * Isp1 * g + 1) / (alpha * Isp1 * g * eps1));
    R2 = ((alpha * Isp2 * g + 1) / (alpha * Isp2 * g * eps2));
    R3 = ((alpha * Isp3 * g + 1) / (alpha * Isp3 * g * eps3));

    % Payload ratios
    lambda1 = (1 - R1 * eps1) / (R1 - 1);
    lambda2 = (1 - R2 * eps2) / (R2 - 1);
    lambda3 = (1 - R3 * eps3) / (R3 - 1);

    % Wet mass of each stage
    m_03 = m_upper * (1 + lambda3) / lambda3;
    m_02 = m_03 * (1 + lambda2) / lambda2;
    m_01 = m_02 * (1 + lambda1) / lambda1;

    % Wet mass of individual stages
    m_3 = m_03 - m_upper;
    m_2 = m_02 - m_03;
    m_1 = m_01 - m_02;

    % Burnout mass of each stage
    m_b1 = m_01 / R1;
    m_b2 = m_02 / R2;
    m_b3 = m_03 / R3;

    % Propellant mass of each stage
    m_p1 = m_01 - m_b1;
    m_p2 = m_02 - m_b2;
    m_p3 = m_03 - m_b3;

    % Structural mass of each stage
    m_s1 = m_b1 - m_02;
    m_s2 = m_b2 - m_03;
    m_s3 = m_b3 - m_upper;
end

function F = optimal(x, Isp1, Isp2, Isp3)
    % x(1): delta u1
    % x(2): delta u2
    % x(3): delta u3
    % x(4): alpha

    % params
    deltaV = 8000; 
    g = 9.806;

    % Calculate epsilons
    eps1 = calc_eps(Isp1, x(1));
    eps2 = calc_eps(Isp2, x(2));
    eps3 = calc_eps(Isp3, x(3));

    % delta ui = Isp * g * ln(Ri)
    % delta ui = Isp * g * ln((alpha * Isp * g + 1) / (alpha * Isp * g * eps))
    F(1) = x(1) - Isp1 * g * log((x(4) * Isp1 * g + 1) / (x(4) * Isp1 * g * eps1));
    F(2) = x(2) - Isp2 * g * log((x(4) * Isp2 * g + 1) / (x(4) * Isp2 * g * eps2));
    F(3) = x(3) - Isp3 * g * log((x(4) * Isp3 * g + 1) / (x(4) * Isp3 * g * eps3));

    % sum of each delta ui is deltaV
    F(4) = deltaV - x(1) - x(2) - x(3);
end

function eps = calc_eps(Isp, delta_u)
    % Calculate epsilon based on Isp
    if (Isp == 460) % hydrogen
        eps = 0.07060945 + 0.1610852 * exp(-0.849934 * (0.001 * delta_u));
    else % hydrocarbon
        eps = 0.0305466 + 0.06892734 * exp(-0.8586846 * (0.001 * delta_u));
    end
end
