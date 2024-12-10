function out = solveCDNozzle(in)
    %%% Required Input Structure %%%
    % in.T0 = 3200;      % stagnation temp (in K)
    % in.P0 = 65;        % stagnation pressure (in atm)
    % in.gamma = 1.2;    % specific heat ratio
    % in.AcAt = 5;       % area ratio of combustion chamber to throat
    % in.AeAt = 15;      % area ratio of exit to throat
    % in.At = 0.01;      % throat area (in m^2), (NOT NEEDED for Q1)
    % in.MW = 13;        % molecular weight (in g/mol)
    % in.Pa = 1;         % ambient pressure (in atm)

    %%% Output structure %%%
    % out.flow : flow type
    % out.P : Pressure inside the nozzle (in atm)
    % out.T : Temperature inside the nozzle (in K)
    % out.M : Mach number inside the nozzle
    % out.u : Axial velocity inside the nozzle (in m/s)
    % out.cF : Thrust coefficient
    % out.Isp : Specific impulse (in s)

    % isoentropic Mach numbers at the exit - subsonic and supersonic solutions
    AxAt = in.AeAt;
    gm = in.gamma;
    areaMach = @(M) areaMach_fcn(M, AxAt, in);

    Me_sub = fsolve(areaMach,0.5,optimoptions('fsolve','Display','none'));   % subsonic solution
    Pe_sub = getPres(Me_sub, in); % exit press at subsonic solution

    Me_sup = fsolve(areaMach,1.5,optimoptions('fsolve','Display','none'));   % supersonic solution
    Pe_sup = getPres(Me_sup, in); % exit press at supersonic solution

    Me_NS = ((1+((gm-1)/2)*Me_sup^2)/(gm*Me_sup^2 - ((gm-1)/2)))^(0.5); % normal shock at exit
    Pe_NS = getPresNS(Me_NS, in); % exit pressure if normal shock at the exit

    Pa = in.Pa;
    % Types of flow
    if Pa > Pe_sub
        out.flow = 'Fully Subsonic';
        [out.P,out.T,out.M,out.u] = backsolveSubsonic(in);
    elseif Pa == Pe_sub
        out.flow = 'Critically Chocked';
        [out.P,out.T,out.M,out.u] = backsolveSubsonic(in);
    elseif Pa < Pe_sub && Pa > Pe_NS
        out.flow = 'Normal Shock';
        [out.P,out.T,out.M,out.u,out.xL_NS,out.Me] = solveForNormalShock(in);
    elseif Pa == Pe_NS
        out.flow = 'Normal Shock at exit';
        [out.P,out.T,out.M,out.u,out.xL_NS,out.Me] = solveForNormalShock(in);
    elseif Pa < Pe_NS && Pa > Pe_sup
        out.flow = 'Over-expanded';
        [out.P,out.T,out.M,out.u] = solveForSupersonicIso(in);
    elseif Pa == Pe_sup
        out.flow = 'Perfectly Expanded';
        [out.P,out.T,out.M,out.u] = solveForSupersonicIso(in);
    elseif Pa < Pe_sup
        out.flow = 'Under-expanded';
        [out.P,out.T,out.M,out.u] = solveForSupersonicIso(in);
    else
        out.flow = 'ERROR';
    end

    % get performance parameters
    [out.cF, out.Isp] = getPerformanceParams(out, in);
end

function [cF, Isp] = getPerformanceParams(out, in)
    gm = in.gamma; P0 = in.P0; Pa = in.Pa; AeAt = in.AeAt; T0 = in.T0;
    MW = in.MW*0.001; % molecular weight in kg/mol
    R = 8.3144; % universal gas constant in J/mol-K
    g = 9.8066; % acceleration due to gravity in m/s2

    % exit conditions
    Pe = out.P(201);

    % thrust coefficient (cF)
    cF = sqrt(((2*gm^2)/(gm-1))*(2/(gm+1))^((gm+1)/(gm-1)) * ...
        (1 - (Pe/P0)^((gm-1)/gm))) + (Pa/P0)*(Pe/Pa - 1)*(AeAt);

    % specific impulse (Isp)
    Isp = (1/g)*sqrt((2*gm*R*T0/(MW*(gm-1)))*(1-(Pe/P0)^((gm-1)/gm))) + ...
          (1/g)*sqrt((1/gm)*(R/MW)*((2/(gm+1))^(-(gm+1)/(gm-1)))) * ...
          (Pa/P0)*(Pe/Pa-1)*sqrt(T0)*(AeAt);
end

function [P,T,M,u] = backsolveSubsonic(in)
    P0 = in.P0; T0 = in.T0; AeAt = in.AeAt; gm = in.gamma; Pa = in.Pa; AcAt = in.AcAt;
    R = 8.3144; % universal gas constant

    Pe = Pa;
    Me = getMachFromPres(Pe, in);

    eq_A_star = @(AxAt) areaMach_fcn(Me, AxAt, in);
    AeA_star = fsolve(eq_A_star,2,optimoptions('fsolve','Display','none')); % area ratio with virtual throat

    % 201 data points from -1 to 1 (x/L)
    P = zeros(201,1); T = zeros(201,1);
    M = zeros(201,1); u = zeros(201,1);

    for i = 1:201
        k = 202 - i;
        xL = -1+((k-1)/100);
        if k >= 101
            AAt = 1+(AeAt-1)*xL;
        else
            AAt = 1-(AcAt-1)*xL;
        end

        AA_star = AAt*AeA_star*(AeAt^(-1));
        eq_mach = @(M_var) areaMach_fcn(M_var, AA_star, in);

        if k == 201
            M_sol = fsolve(eq_mach,0.2,optimoptions('fsolve','Display','none')); % subsonic solution
        else
            M_sol = fsolve(eq_mach,M(k+1),optimoptions('fsolve','Display','none')); % previous solution as guess
        end

        P(k) = P0*(1+((gm-1)/2)*M_sol^2)^(-gm/(gm-1));
        T(k) = T0*(1+((gm-1)/2)*M_sol^2);
        M(k) = M_sol;
        u(k) = M(k)*sqrt(gm*R*T(k));
    end
end

function [P,T,M,u] = solveForSupersonicIso(in)
    % converging section
    [P,T,M,u] = solveConvergingIsoentropic(in);

    % diverging section
    [P,T,M,u] = solveDivergingUpstream(P,T,M,u,1,in);
end

function [P,T,M,u,xL_NS,Me] = solveForNormalShock(in)
    gm = in.gamma; Pa = in.Pa; P0 = in.P0; AeAt = in.AeAt;
    Pe = Pa;

    pres_area_ratio = (Pe/P0)*(AeAt);

    % equation for Me (mach number at exit)
    eq_Me = @(Me) (Pe/P0 - (1+((gm-1)/2)*Me^2)^(-gm/(gm-1))*(1/AeAt));
    Me = fsolve(eq_Me,0.5,optimoptions('fsolve','Display','none'));

    % ratio of stagnation pressure to exit pressure
    P0e_Pe_ratio = (1+((gm-1)/2)*Me^2)^(gm/(gm-1));

    % ratio of downstream upstream stagnation pressure
    P0d_P0u_ratio = P0e_Pe_ratio*(Pe/P0);

    % solving for mach number at normal shock
    eq_Ms = @(Ms) -P0d_P0u_ratio + ...
        (1+((gm-1)/2)*((1+((gm-1)/2)*Ms^2)/(gm*Ms^2-((gm-1)/2))))^(gm/(gm-1)) * ...
        (1+((2*gm)/(gm+1))*(Ms^2-1)) * ...
        (1+((gm-1)/2)*Ms^2)^(-gm/(gm-1));
    Ms = fsolve(eq_Ms,2,optimoptions('fsolve','Display','none'));

    % area ratio from mach number
    eq_NS_area_ratio = @(AsAt) areaMach_fcn(Ms, AsAt, in);
    AsAt = fsolve(eq_NS_area_ratio,2,optimoptions('fsolve','Display','none'));

    % location of the normal shock (diverging section)
    xL_NS = (AsAt - 1)/(AeAt - 1);

    % downstream mach number
    Md = ((1+((gm-1)/2)*Ms^2)/(gm*Ms^2 - ((gm-1)/2)))^(1/2);

    % pressure, temperature, mach number, velocity distribution
    [P,T,M,u] = solveConvergingIsoentropic(in);

    % upstream supersonic isoentropic (before normal shock)
    [P,T,M,u] = solveDivergingUpstream(P,T,M,u,xL_NS,in);

    % downstream subsonic isoentropic (after normal shock)
    [P,T,M,u] = solveDivergingDownstream(P,T,M,u,xL_NS,Md,AsAt,in);
end

function [P,T,M,u] = solveDivergingDownstream(P,T,M,u,xL_NS,Md,AsAt,in)
    P0 = in.P0; T0 = in.T0; AeAt = in.AeAt; gm = in.gamma;
    R = 8.3144; % universal gas constant

    eq_Ad_star = @(AxAt) areaMach_fcn(Md,AxAt,in);
    AsAd_star = fsolve(eq_Ad_star,2,optimoptions('fsolve','Display','none')); % area ratio with virtual throat

    P0d = P0*(AsAd_star/AsAt);

    for k = fix(100*xL_NS)+102:201
        xL = -1+((k-1)/100);
        AAt = 1+(AeAt-1)*xL;
        AAd_star = AAt*AsAd_star*(AsAt^(-1));
        eq_mach = @(M_var) areaMach_fcn(M_var, AAd_star, in);

        if k == fix(100*xL_NS)+102
            M_sol = fsolve(eq_mach,0.2,optimoptions('fsolve','Display','none')); % subsonic solution
        else
            M_sol = fsolve(eq_mach,M(k-1),optimoptions('fsolve','Display','none')); % previous solution as guess
        end

        P(k) = P0d*(1+((gm-1)/2)*M_sol^2)^(-gm/(gm-1));
        T(k) = T0*(1+((gm-1)/2)*M_sol^2);
        M(k) = M_sol;
        u(k) = M(k)*sqrt(gm*R*T(k));
    end
end

function [P,T,M,u] = solveConvergingIsoentropic(in)
    % 201 data points from -1 to 1 (x/L)
    P = zeros(201,1); T = zeros(201,1);
    M = zeros(201,1); u = zeros(201,1);

    P0 = in.P0; T0 = in.T0; AcAt = in.AcAt; gm = in.gamma;
    R = 8.3144; % universal gas constant

    for k = 1:101
        xL = -1+((k-1)/100);
        AAt = 1-(AcAt-1)*xL;
        eq_mach = @(M_var) areaMach_fcn(M_var,AAt,in);

        if k == 1
            M_sol = fsolve(eq_mach,0.2,optimoptions('fsolve','Display','none')); % subsonic solution
        else
            M_sol = fsolve(eq_mach,M(k-1),optimoptions('fsolve','Display','none'));
        end

        P(k) = P0*(1+((gm-1)/2)*M_sol^2)^(-gm/(gm-1));
        T(k) = T0*(1+((gm-1)/2)*M_sol^2);
        M(k) = M_sol;
        u(k) = M(k)*sqrt(gm*R*T(k));
    end
end

function [P,T,M,u] = solveDivergingUpstream(P,T,M,u,xL_NS,in)
    P0 = in.P0; T0 = in.T0; AeAt = in.AeAt; gm = in.gamma;
    R = 8.3144; % universal gas constant

    for k = 102:fix(100*xL_NS)+101
        xL = -1+((k-1)/100);
        AAt = 1+(AeAt-1)*xL;
        eq_mach = @(M_var) areaMach_fcn(M_var, AAt, in);

        if k == 102
            M_sol = fsolve(eq_mach,1.2,optimoptions('fsolve','Display','none')); % supersonic solution
        else
            M_sol = fsolve(eq_mach,M(k-1),optimoptions('fsolve','Display','none'));
        end

        P(k) = P0*(1+((gm-1)/2)*M_sol^2)^(-gm/(gm-1));
        T(k) = T0*(1+((gm-1)/2)*M_sol^2);
        M(k) = M_sol;
        u(k) = M(k)*sqrt(gm*R*T(k));
    end
end

% Mach number as a function of pressure for isoentropic flow
function M_sol = getMachFromPres(Pe,in)
    gm = in.gamma; P0 = in.P0;
    eq = @(M) (Pe/P0 - (1+((gm-1)/2)*M^2)^(-gm/(gm-1)));
    M_sol = fsolve(eq,0.5,optimoptions('fsolve','Display','none')); % subsonic solution
end

% Pressure as a function of Mach number
function F = getPres(M, in)
    P0 = in.P0; gm = in.gamma;
    F = P0*(1+((gm-1)/2)*M^2)^(-gm/(gm-1));
end

% Pressure as a function of Mach number after Normal Shock
function F = getPresNS(M,in)
    P0 = in.P0; gm = in.gamma; AeAt = in.AeAt;

    eq_Ad_star = @(AxAt) areaMach_fcn(M,AxAt,in);
    AdAd_star = fsolve(eq_Ad_star,2,optimoptions('fsolve','Display','none')); % area ratio with virtual throat

    P0d = P0*(AdAd_star/AeAt);
    F = P0d*(1+((gm-1)/2)*M^2)^(-gm/(gm-1));
end

% Temperature as a function of Mach number
function F = getTemp(M,in)
    T0 = in.T0; gm = in.gamma;
    F = T0*(1+((gm-1)/2)*M^2);
end

% Area-Mach relations for choked flow (Mt = 1 assumption)
function F = areaMach_fcn(M,AxAt,in)
    gm = in.gamma;
    F = AxAt^2 - (1/(M^2))*((2/(1+gm))*(1+((gm-1)/2)*M^2))^((gm+1)/(gm-1));
end

%% Test Code

clear all; close all; clc;

% Input parameters structure (default values)
in.T0 = 3200;    % Stagnation temp in K
in.P0 = 65;      % Stagnation pressure in atm
in.gamma = 1.2;  % Specific heat ratio
in.AcAt = 5;     % Area ratio of combustion chamber to throat
in.AeAt = 15;    % Area ratio of exit to throat
in.At = 0.01;    % Throat area in m^2 (not needed for Q1)
in.MW = 13;      % Molecular weight of propellant product in g/mol
in.Pa = 1;       % Ambient pressure in atm (initially)

% Problem 1: Numerical Examples at select conditions
out = cell(5,1);

in.Pa = 1;        % Sea-level
out{1} = solveCDNozzle(in);

in.Pa = 0.485;    % ~6 km altitude
out{2} = solveCDNozzle(in);

in.Pa = 0.12;     % ~15 km altitude
out{3} = solveCDNozzle(in);

in.Pa = 0.000750062; % ~50 km altitude (≈76 Pa)
out{4} = solveCDNozzle(in);

in.Pa = 15;       % Example back pressure for Normal Shock
out{5} = solveCDNozzle(in);

% Set up plotting
x_L = linspace(-1, 1, 201);
figure('Position',[10 10 1500 1500]);

% Plotting each parameter for all cases
% Pressure
subplot(2,2,1);
for i = 1:5
    plot(x_L, out{i}.P / in.P0, 'LineWidth', 2); hold on;
end
xlabel('x / L');
ylabel('p(x) / p0');
title('p(x) / p0 vs. x / L');

% Temperature
subplot(2,2,2);
for i = 1:5
    plot(x_L, out{i}.T / in.T0, 'LineWidth', 2); hold on;
end
xlabel('x / L');
ylabel('T(x) / T0');
title('T(x) / T0 vs. x / L');

% Mach Number
subplot(2,2,3);
for i = 1:5
    plot(x_L, out{i}.M, 'LineWidth', 2); hold on;
end
xlabel('x / L');
ylabel('M(x)');
title('M(x) vs. x / L');

% Velocity
subplot(2,2,4);
for i = 1:5
    plot(x_L, out{i}.u, 'LineWidth', 2); hold on;
end
xlabel('x / L');
ylabel('u(x)');
title('u(x) vs. x / L');

% Add legends
subplot(2,2,1)
legend('P_a = 1 atm (over-expanded)', 'P_a = 0.485 atm (over-expanded)', ...
       'P_a = 0.12 atm (under-expanded)', 'P_a = 76 Pa (under-expanded)', ...
       'P_a = 15 atm (Normal Shock)', 'Location', 'best')
grid on

subplot(2,2,2)
legend('P_a = 1 atm (over-expanded)', 'P_a = 0.485 atm (over-expanded)', ...
       'P_a = 0.12 atm (under-expanded)', 'P_a = 76 Pa (under-expanded)', ...
       'P_a = 15 atm (Normal Shock)', 'Location', 'northwest')
grid on

subplot(2,2,3)
legend('P_a = 1 atm (over-expanded)', 'P_a = 0.485 atm (over-expanded)', ...
       'P_a = 0.12 atm (under-expanded)', 'P_a = 76 Pa (under-expanded)', ...
       'P_a = 15 atm (Normal Shock)', 'Location', 'best')
grid on

subplot(2,2,4)
legend('P_a = 1 atm (over-expanded)', 'P_a = 0.485 atm (over-expanded)', ...
       'P_a = 0.12 atm (under-expanded)', 'P_a = 76 Pa (under-expanded)', ...
       'P_a = 15 atm (Normal Shock)', 'Location', 'best')
grid on