function [Me, Te, Pe, ue, Ae_At, F, Isp, xw, yw, trunc_index, Prat, Xmesh, noz_len_1m ,noz_ex_rad_1m]=calcEngine(in, optimal_mode, trunc_length, scale) % chamber pressure & temperature are specified as input parameters. These are constrained by current technology limits.
    % We assume the engine is optimized at either sea level and vacuum, and ignore off-point of operation situations like during
    % the ascend phase of first stage, transition between sea level and vacuum while first stage is only optimize for sea level. 
    % Assuming everything downstream of combustion chamber 
    % (converging section, throat, diverging section) are isentropic (no shock presence, minimal heat transfer b/w exhaust and wall)
    gamma = in.gamma;
    T0 = in.T0;
    P0 = in.P0; % Chamber condition upon combustion of propellants, known values from the beginning
    At = in.At;
    R_specific = in.R_specific;
    %-------------Analysis------------%

    if strcmp(optimal_mode, 'sea level')
        Pe = 700000; % Sea level exit pressure
    elseif strcmp(optimal_mode, 'vacuum')
        Pe = 200000; % vacuum exit pressure, since it's impossible to really expand it to zero pressure, we take 1000 pa as required exit pressure      
    else
        fprintf('error');
    end


    Me = Mach_from_expansion(Pe, P0, gamma); % Exit Mach number
    % Calculate Full Nozzle and then truncate it
    trunc_length = trunc_length * scale; % Bring it to the same scale
    [xw,yw,xcl,Mcl,Prat, X, Y] = MinLenNozDes(scale, Me,in.gamma,80,0.1, 0);
    Xmesh = X;
    % Find the x index where elements exceed the truncation length for the last time
    
    for i = 1:length(xw)
        if xw(i) >= xw(end) - trunc_length
            trunc_index = i;
            break
        end
    end
    for i = 1:length(X(1, :))
        if X(1, i) >= X(1, end) - trunc_length
            mesh_index = i;
            break
        end
    end

    % Check if truncation_index is empty, which can occur if all values exceed truncation_length
    Pe = mean(P0*Prat(:, mesh_index),'omitnan');
    Me = Mach_from_expansion(Pe, P0, gamma); % Now we compute the exit mach number given the truncated exit pressure


    clf
    fig = figure();
    pcolor(X(:, 1:mesh_index),Y(:, 1:mesh_index),Prat(:, 1:mesh_index))
    hold on
    pcolor(X(:, 1:mesh_index),-Y(:, 1:mesh_index),Prat(:, 1:mesh_index))
    plot(xw(1:trunc_index),yw(1:trunc_index),'k-','LineWidth',2) 
    plot(xw(1:trunc_index),-yw(1:trunc_index),'k-','LineWidth',2)
    shading interp
    formatplot('$p / p_o$')
    uiwait(fig)

    Te_T0 = temperature_ratio(Me, gamma); % Exit/Chamber pressure ratio
    Te = Te_T0 * T0; % Exit pressure
    ue = velocity(Me, gamma, R_specific, Te); % Exit velocity
    Ae_At = A_At_from_Mach(Me, gamma); % Nozzle exit-throat area ratio

    rho0 = P0/(R_specific*T0); % Stagnation density
    m_dot = (density_ratio(1, gamma) * rho0) * At * sqrt(gamma * R_specific * T0 * temperature_ratio(1, gamma)); % Calculate mass flow rate at using throat conditions
    if strcmp(optimal_mode, 'sea level')
        Pa = 101325; % Sea level ambient pressure
    elseif strcmp(optimal_mode, 'vacuum')
        Pa = 0; % vacuum ambient pressure, since it's impossible to really expand it to zero pressure, we take 1000 pa as required exit pressure      
    else
        fprintf('error')
    end
    F = m_dot * ue + (Pe - Pa) * Ae_At * At; % Calculate thrust of the nozzle flow
    Isp = F/(m_dot * 9.806);
    noz_len_1m = xcl(end)/scale; % Nozzle length if the throat is 1m in radius
    noz_ex_rad_1m = yw(end)/scale; % Nozzle exit radius  if the throat is 1m in radius
end 


function Ma = area_ratio_to_mach(A_At, gamma)
    % Solves for Mach number (Ma) given area ratio (A/A*) and specific heat ratio (gamma).

    % Define the equation to solve
    func = @(M) (1 ./ M) .* ((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M.^2)).^((gamma + 1) / (2 * (gamma - 1))) - A_At;

    % Initial guess for Mach number
    initial_guess = 2;  % Starting guess for supersonic flow

    % Solve the equation using fsolve
    options = optimoptions('fsolve', 'Display', 'off');  % Turn off display for solver
    Ma = fsolve(func, initial_guess, options);

end

function pr = pressure_ratio(M, gamma)
    % Calculates the pressure ratio across a shock wave given Mach number (M) and specific heat ratio (gamma)
    pr = (1 + (gamma - 1) / 2 * M.^2).^(-gamma / (gamma - 1));
end

function dr = density_ratio(M, gamma)
    % Calculates the density ratio across a shock wave given Mach number (M) and specific heat ratio (gamma)
    dr = (1 + (gamma - 1) / 2 * M.^2).^(-1 / (gamma - 1));
end

function tr = temperature_ratio(M, gamma)
    % Calculates the temperature ratio across a shock wave given Mach number (M) and specific heat ratio (gamma)
    tr = (1 + (gamma - 1) / 2 * M.^2).^(-1);
end

function v = velocity(M, gamma, Rspecific, T)
    % Calculates the flow velocity given Mach number (M), specific heat ratio (gamma), specific gas constant (Rspecific), and static temperature (T)
    v = M .* sqrt(gamma * Rspecific * T);
end

function Ma = Mach_from_expansion(P, P0, gamma)
    % Calculate the Mach number from given pressures and specific heat ratio
    pressureRatio = P / P0;
    exponent = -(gamma - 1) / gamma;
    Ma = sqrt((2 / (gamma - 1)) * (pressureRatio^exponent - 1));
end

function A_At = A_At_from_Mach(M, gamma)
    % Calculate the area ratio A/A* for a given Mach number M and specific heat ratio gamma
    term1 = 1 / M;
    term2 = ((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M^2)) ^ ((gamma + 1) / (2 * (gamma - 1)));
    A_At = term1 * term2;  % Area of shock over area of throat
end

function formatplot(clab)
    xlabel('Distance from throat, mm', 'interpreter', 'latex', 'fontsize',24,'fontname','times')
    ylabel('Distance from axis, mm',   'interpreter', 'latex', 'fontsize',24,'fontname','times')
    
    set(gca,'fontsize',24,'fontname','times')
    set(gca,'linewidth',1.5,'box','off','ticklength',[.01 0])
    set(gca,'tickdir','out')
    
    axis equal
    
    if ~isempty(clab)
        h = colorbar;
        ylabel(h, clab,'fontsize',24,'fontname','times','interpreter', 'latex')
    end
    
    set(gcf,'position',[100 100 900 600],'color',[1 1 1],'paperPositionMode','auto')
    set(gca, 'position', [0.06,0.12,0.88, 0.85])
        
    end