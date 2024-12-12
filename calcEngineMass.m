function mass = calcEngineMass(in, rho_wall, R_c, L_c, sigma_h, nlines, thi, oflag)
    yt = in.At / (2 * in.throat_width);
    disp(yt)
    disp(in.Me)
    disp(in.gamma)
    disp(nlines)
    disp(thi)
    disp(oflag)
    [xw,yw,xcl,Mcl] = MinLenNozDes(yt, in.Me,in.gamma,nlines,thi,oflag);
    % [xw,yw,xcl,Mcl] = MinLenNozDes(0.3536, 2.4475, 1.1443, 20, 0.1, 1)
    p_c = in.P0;

    wall_thickness = p_c * R_c / sigma_h;
    disp(wall_thickness)
    nozzle_material_volume = 0; % volume occupied by the wall material of the chamber (volume of the shell excluding its interior)
    for i = 1:length(xw)-1
        % TO DO: Add mach number & pressure code to decrease thickness as
        % it goes
        nozzle_width = (yw(i) + yw(i+1))/2;
        nozzle_material_volume = nozzle_material_volume + wall_thickness * (xw(i+1) - xw(i)) * nozzle_width * 4;
    end
    mass_nozzle = nozzle_material_volume * rho_wall;
    mass_chamber = 2 * pi * R_c * L_c * wall_thickness * rho_wall; % Assume cylindrical chamber, with chamber radius and height known. 
    % Assume thin wall (thus ignore variation of inner and outer diameter)
    mass = mass_chamber + mass_nozzle;
end