function F = SolarEnergyBalances(x,y,z,wLevel,sqArea)
    %EVERYTHING AS HEAT FLUX CURRENTLY
    
    solar_flux = y;   %fnc takes solar flux as input
    Ta = z;        %fnc takes ambient air temp as input
    Tsky = 0.0552*Ta^(1.5);  %sky temperature
    
    %Radiation
    SBConst = 5.67*10^(-8);  %Stefan Boltzmann constant
    g_emis = 0.95;              %glass emissivity
    w_emis = 0.96;              %water emissivity
    
    %Transmissivity
    g_transm = 0.92;       %transmittance of glass
    w_transm = 0.9;         %transmittance of water
    
    %Absorbtivity
    g_absorb = 0.016;         %absorbance of glass
    w_absorb = 0.0*g_transm;
    b_absorb = 0.84*g_transm*w_transm;
    b_absorb = 0.7*g_transm*w_transm;
    
    %Geometry
    wallHeight = 0.04;      %height of solar still walls
    wallThick = 0.01;       %thickness of solar still walls
    coverAngle = 20;        %angle the glass makes with respect to horizon
    A = sqArea;                  %water and basin surface area
    Ag = 2*(sqrt(A)*(sqrt(A)*0.5)/cosd(coverAngle));     %glass surface area
    Ains = sqrt(A) + 4*wallHeight*sqrt(A);
    
  
    %Hc_ga = 10;             %convection between glass and air
    Hc_ga = 20;
    Hc_ia = 2;
    
    Kw = 0.6;               %conduction coef for water
    Lw = wLevel;            %height of water
    %Ki = 0.06;
    Ki = 0.04;              %conduction coef of insulation
    Li = 0.02;              %height of insulation
    %Li = 0.01;

    
    Pw = exp(25.317-5144/(273.15+x(2)));    %partial pressure of vapor at water
    Pg = exp(25.317-5144/(273.15+x(1)));    %partial pressure of vapor at glass

    Hc_gw = 0.884*(x(2)-x(1) + x(2)*(Pw-Pg)/(268.9*1000 - Pw))^(1/3);   %convection coef between water and glass
    He_gw = (16.273*10^(-3))*Hc_gw*(Pw-Pg)/(x(2)-x(1));                 %evaporative coef between water and glass

    %glass energy balance
    F(1)=(Ag*solar_flux*g_absorb + A*Hc_gw*(x(2)-x(1)) + A*He_gw*(x(2)-x(1)) + A*SBConst*w_emis*(((x(2)+273.15)^4)-((x(1)+273.15)^4)) - Ag*Hc_ga*(x(1)-Ta) - Ag*SBConst*g_emis*(((x(1)+273.15)^4)-(Tsky+273.15)^4));
    
    %water energy balance
    F(2)=(A*solar_flux*w_absorb + A*(Kw/Lw)*(x(3)-x(2)) - A*Hc_gw*(x(2)-x(1)) - A*He_gw*(x(2)-x(1)) - A*SBConst*w_emis*(((x(2)+273.15)^4)-((x(1)+273.15)^4)));
    
    %basin energy balance
    %F(3)=(A*solar_flux*b_absorb - Ains*(Ki/Li)*(((x(3)+x(2)+x(1))/3)-Ta) - A*(Kw/Lw)*(x(3)-x(2)));
    %F(3)=(A*solar_flux*b_absorb - Ains*(Ki/Li)*(((x(3)+x(2)+x(1))/3)-x(4)) - A*(Kw/Lw)*(x(3)-x(2)));
    %F(3)=(A*solar_flux*b_absorb - Ains*(Ki/Li)*(x(3)-x(4)) - A*(Kw/Lw)*(x(3)-x(2)));
    F(3)=(A*solar_flux*b_absorb - Ains*(Ki/Li)*(x(3)-((Ta+x(3))/2)) - A*(Kw/Lw)*(x(3)-x(2)));

    
    
  
    