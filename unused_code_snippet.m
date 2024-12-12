% After all engine numbers are obtained, assign Isp to each stage based on engine selection
if strcmp(fuel_type_stg1, 'hydrogen')
    Isp1 = Isp_HY;
else
    Isp1 = Isp_HC;
end

if strcmp(fuel_type_stg2, 'hydrogen')
    Isp2 = Isp_HY;
else
    Isp2 = Isp_HC;
end

if strcmp(fuel_type_stg3, 'hydrogen')
    Isp3 = Isp_HY;
else
    Isp3 = Isp_HC;
end


