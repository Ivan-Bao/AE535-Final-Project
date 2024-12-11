function mass = calcEngineMass(in, rho_wall, R_c, L_c, sigma_h, nlines, thi, oflag)
    yt = in_At / (2 * in.throat_width);
    [xw,yw,xcl,Mcl] = MinLenNozDes(yt,in.Me,in.gamma,nlines,thi,oflag,pflag);
    p_c = in.P0;

    wall_thickness = p_c * R_c / sigma_h;
    nozzle_material_volume = 0; % volume occupied by the wall material of the chamber (volume of the shell excluding its interior)
    for i = 1:length(xw)-1
        nozzle_width = (yw(i) + yw(i+1))/2;
        nozzle_material_volume = nozzle_material_volume + wall_thickness * (xw(i+1) - xw(i)) * nozzle_width * 4;
    end
    mass_nozzle = nozzle_material_volume * rho_wall;
    mass_chamber = 2 * pi * R_c * L_c * wall_thickness * rho_wall;
    mass = mass_chamber + mass_nozzle;
end