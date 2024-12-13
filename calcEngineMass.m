function [mass, Prat] = calcEngineMass(in, trunc_index, rho_wall, R_c, L_c, sigma_h, xw, yw, scale, Prat, Xmesh)
    p_c = in.P0;

    wall_thickness = p_c * R_c / (sigma_h * scale);
    mass_chamber = 2 * pi * R_c * L_c * wall_thickness * rho_wall; % Assume cylindrical chamber, with chamber radius and height known. 
    nozzle_material_volume = 0; % volume occupied by the wall material of the chamber (volume of the shell excluding its interior)
    for i = 1:trunc_index-1

        [~, mesh_index] = min(abs(Xmesh(1,:) - xw(i)));
        new_wall_thickness = p_c * Prat(round(0.5*length(Prat(:, 1))), mesh_index)*  R_c/ (2 * sigma_h);
        if not(isnan(new_wall_thickness))
            wall_thickness = new_wall_thickness; % if we are not having nan data point of pressure expansion
        end
        nozzle_width = (yw(i) + yw(i+1))/(2*scale);
        nozzle_material_volume = nozzle_material_volume + wall_thickness * (xw(i+1) - xw(i)) * nozzle_width * 4 / scale;
    end
    mass_nozzle = nozzle_material_volume * rho_wall;
    
    % Assume thin wall (thus ignore variation of inner and outer diameter)
    mass = mass_chamber + mass_nozzle;
end