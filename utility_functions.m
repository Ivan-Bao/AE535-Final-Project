
classdef utility_functions
    methods(Static)
        function epsilon = structural_coefficient(fuel_type, delta_v)
            % Calculate the structural coefficient based on the fuel type and delta-v.
            %
            % Parameters:
            % - fuel_type: either 'hydrogen' or 'hydrocarbon'
            % - delta_v: change in velocity in m/s
            %
            % Returns:
            % - epsilon: the structural coefficient
        
            % Placeholder coefficients for exponential fit: a, b, c
            % These should be replaced with actual values from your curve fitting.
            if strcmp(fuel_type, 'hydrogen')
                a = 1.0; % Example coefficient a for hydrogen
                b = 0.0005; % Example coefficient b for hydrogen
                c = 0.1; % Example coefficient c for hydrogen
            elseif strcmp(fuel_type, 'hydrocarbon')
                a = 1.2; % Example coefficient a for hydrocarbon
                b = 0.0003; % Example coefficient b for hydrocarbon
                c = 0.2; % Example coefficient c for hydrocarbon
            else
                error('Invalid fuel type. Choose ''hydrogen'' or ''hydrocarbon''');
            end
        
            % Exponential function to model the structural coefficient
            epsilon = a * exp(-b * delta_v) + c;
        end
    end
end
