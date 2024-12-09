
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
                epsilon = 0.07014715 + 0.1648617* exp ( -0.8647945*(0.001*delta_v));
            elseif strcmp(fuel_type, 'hydrocarbon')
                epsilon = 0.0305466 + 0.06892734* exp ( -0.8586846*(0.001*delta_v))

            else
                error('Invalid fuel type. Choose ''hydrogen'' or ''hydrocarbon''');
            end

        end
    end
end
