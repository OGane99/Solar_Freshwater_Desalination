clear all
close all
clc

timestepLength = 0.25; %steplength of one is 1 hour long, 0.5 is 30 min etc..
amountSteps = 24/timestepLength;    % total amount of time steps required 

Daily_Output = zeros([amountSteps,1]);
Time_Matrix = zeros([amountSteps,1]);
Efficiency_Matrix = zeros([amountSteps,1]);
Irradiation_Matrix = zeros([amountSteps,1]);
%Savethis_Matrix = zeros([amountSteps,1]);


for i = 1:amountSteps
    Time_Matrix([i,1])= i*timestepLength-1;
    Time_Matrix([1,1])=0;
end


%initialize variables that cannot be included in energy balance
Gmax = 460;  %max solar irradiance in worst conditions
Gmax = 890;  %max solar irradiance in worst conditions
%Gmax = (890+460)/2;

latentHeatVaporization = 2.25*10^6;
T_amb0 = 20;                    %initial ambient temperature
startingWaterLevel = 0.02;      %initial water level
A = 1.5;                  %water and basin surface area (only for water level changes in this code)


for i = 1:amountSteps
    if i == 1
        wLevel = startingWaterLevel;
    else
        wLevel = wLevel - (Litres_hourly_clean_water/A)/1000;
    end
    
    %Solar irradiation as a function of time (hrs)
    G = (Gmax/2)*(sin(pi*(i*timestepLength-1)/11.51))+abs((Gmax/2)*(sin(pi*(i*timestepLength-1)/11.51)));
    %G = (Gmax/2)*(sin(pi*(i*timestepLength-1)/9.83))+abs((Gmax/2)*(sin(pi*(i*timestepLength-1)/9.83)));
    %G = (Gmax/2)*(sin(pi*(i*timestepLength-1)/8.1588))+abs((Gmax/2)*(sin(pi*(i*timestepLength-1)/8.1588)));
    Irradiation_Matrix([i,1]) = G*timestepLength;
    
    Tamb = T_amb0 + T_amb0*0.5*sin(pi*i*timestepLength/24);

    %solve the system of nonlinear equations
    %with initial guess of Tw, Tg, Tb
    fun = @(x) SolarEnergyBalances(x, G, Tamb, wLevel,A); 
    x0 = [20,21,22];    
    x = fsolve(fun,x0); 

    %from solution, use Tw to find evaporative heat transfer, and
    %from that, hourly water output
    Pw = exp(25.317-5144/(273.15+x(2)));
    Pg = exp(25.317-5144/(273.15+x(1)));
    Hc_gw = 0.884*(x(2)-x(1) + x(2)*(Pw-Pg)/(268.9*1000 - Pw))^(1/3);
    He_gw = (16.273*10^(-3))*Hc_gw*(Pw-Pg)/(x(2)-x(1));
    Litres_hourly_clean_water = A*He_gw*(x(2)-x(1))*(3600*timestepLength)/latentHeatVaporization;

    Daily_Output([i,1])=Litres_hourly_clean_water;
    Efficiency_Matrix([i,1]) = A*He_gw*(x(2)-x(1))/G;
end

%Daily_Output = zeros([amountSteps,1]);
Daily_Output([1,1])=0;

total_G_KWH = sum(Irradiation_Matrix)
total_daily_water_output = sum(Daily_Output)

figure(1);
fprintf('\n\n Solar still daily water output = %1.2f L \n\n',total_daily_water_output);
plot(Time_Matrix,Daily_Output/timestepLength,'b--o','LineWidth',0.7);
hold on;
%plot(Time_Matrix,Daily_Output1/timestepLength,'r--o','LineWidth',0.7)
title('Potable Water Output Over Time for Single Basin Solar Still (Absorptivity = 0.7)')
xlabel('Time since sunrise [hrs]')
ylabel('Potable water output [L/hr]')
xlim([0 12.2])
ylim([-0.1 1.1])
%txt = ['Net Output = ',num2str(round(total_daily_water_output,2)), ' L'];
%text(8,0.9,txt,'FontSize',11);
%legend('Best Net Output = 6.65 L','Worst Net Output = 4.14 L')

figure(2);
plot(Time_Matrix,Efficiency_Matrix,'b--o','LineWidth',0.7);
title('Efficiency of Solar Still in Best and Worst Conditions (Absorptivity = 0.7)')
xlim([0 11.25])
ylim([-0.1 1])
xlabel('Time since sunrise [hrs]')
ylabel('System Efficiency [%]')
hold on;


%%

%clear all
%close all
%clc

timestepLength = 0.25; %steplength of one is 1 hour long, 0.5 is 30 min etc..
amountSteps = 24/timestepLength;    % total amount of time steps required 

Daily_Output1 = zeros([amountSteps,1]);
Time_Matrix = zeros([amountSteps,1]);
Efficiency_Matrix = zeros([amountSteps,1]);
Irradiation_Matrix = zeros([amountSteps,1]);
%Savethis_Matrix = zeros([amountSteps,1]);


for i = 1:amountSteps
    Time_Matrix([i,1])= i*timestepLength-1;
    Time_Matrix([1,1])=0;
end


%initialize variables that cannot be included in energy balance
Gmax = 460;  %max solar irradiance in worst conditions
%Gmax = 890;  %max solar irradiance in worst conditions
%Gmax = (890+460)/2;

latentHeatVaporization = 2.25*10^6;
T_amb0 = 20;                    %initial ambient temperature
startingWaterLevel = 0.02;      %initial water level
A = 1.5;                  %water and basin surface area (only for water level changes in this code)


for i = 1:amountSteps
    if i == 1
        wLevel = startingWaterLevel;
    else
        wLevel = wLevel - (Litres_hourly_clean_water/A)/1000;
    end
    
    %Solar irradiation as a function of time (hrs)
    G = (Gmax/2)*(sin(pi*(i*timestepLength-1)/8.1588))+abs((Gmax/2)*(sin(pi*(i*timestepLength-1)/8.1588)));
    %G = (Gmax/2)*(sin(pi*(i*timestepLength-1)/9.83))+abs((Gmax/2)*(sin(pi*(i*timestepLength-1)/9.83)));
    Irradiation_Matrix([i,1]) = G*timestepLength;
    
    Tamb = T_amb0 + T_amb0*sin(pi*i*timestepLength/24);

    %solve the system of nonlinear equations
    %with initial guess of Tw, Tg, Tb
    fun = @(x) SolarEnergyBalances(x, G, Tamb, wLevel,A); 
    x0 = [20,21,22,23];    
    x = fsolve(fun,x0); 

    %from solution, use Tw to find evaporative heat transfer, and
    %from that, hourly water output
    Pw = exp(25.317-5144/(273.15+x(2)));
    Pg = exp(25.317-5144/(273.15+x(1)));
    Hc_gw = 0.884*(x(2)-x(1) + x(2)*(Pw-Pg)/(268.9*1000 - Pw))^(1/3);
    He_gw = (16.273*10^(-3))*Hc_gw*(Pw-Pg)/(x(2)-x(1));
    Litres_hourly_clean_water = A*He_gw*(x(2)-x(1))*(3600*timestepLength)/latentHeatVaporization;

    Daily_Output1([i,1])=Litres_hourly_clean_water;
    Efficiency_Matrix([i,1]) = A*He_gw*(x(2)-x(1))/G;
end

%Daily_Output = zeros([amountSteps,1]);
Daily_Output1([1,1])=0;

total_G_KWH = sum(Irradiation_Matrix)
total_daily_water_output1 = sum(Daily_Output1)


for i = 1:amountSteps
    %Savethis_Matrix([i,1]) = Daily_Output([i,1]);
    %Savethis_Matrix([i,1]);
end

%fprintf('\n\n Solar still daily water output = %1.2f L \n\n',total_daily_water_output1);
figure(1);
plot(Time_Matrix,Daily_Output1/timestepLength,'r--o','LineWidth',0.7)
legend('Best Net Output = 6.74 L','Worst Net Output = 4.15 L')
%legend('Best Net Output = 5.14 L','Worst Net Output = 2.93 L')
%title('Potable Water Output Over Time for Single Basin Solar Still')
%xlabel('Time since sunrise [hrs]')
%ylabel('Potable water output [L/hr]')
%xlim([0 12.2])
%ylim([-0.1 1])
%legend('Best Net Output = 6.29L','Worst Net Output = 4.13L')

figure(2);
plot(Time_Matrix,Efficiency_Matrix,'r--o','LineWidth',0.7)
legend('Efficiency of best output','Efficiency of worst output')

%%

%clear all
%close all
%clc

timestepLength = 0.25; %steplength of one is 1 hour long, 0.5 is 30 min etc..
amountSteps = 24/timestepLength;    % total amount of time steps required 

Daily_Output2 = zeros([amountSteps,1]);
Time_Matrix = zeros([amountSteps,1]);
Efficiency_Matrix = zeros([amountSteps,1]);
Irradiation_Matrix = zeros([amountSteps,1]);
%Savethis_Matrix = zeros([amountSteps,1]);


for i = 1:amountSteps
    Time_Matrix([i,1])= i*timestepLength-1;
    Time_Matrix([1,1])=0;
end


%initialize variables that cannot be included in energy balance
Gmax = 460;  %max solar irradiance in worst conditions
%Gmax = 890;  %max solar irradiance in worst conditions
Gmax = (890+460)/2;

latentHeatVaporization = 2.25*10^6;
T_amb0 = 25;                    %initial ambient temperature
startingWaterLevel = 0.02;      %initial water level
A = 1.5;                  %water and basin surface area (only for water level changes in this code)


for i = 1:amountSteps
    if i == 1
        wLevel = startingWaterLevel;
    else
        wLevel = wLevel - (Litres_hourly_clean_water/A)/1000;
    end
    
    %Solar irradiation as a function of time (hrs)
    G = (Gmax/2)*(sin(pi*(i*timestepLength-1)/9.83))+abs((Gmax/2)*(sin(pi*(i*timestepLength-1)/9.83)));
    Irradiation_Matrix([i,1]) = G*timestepLength;
    
    Tamb = T_amb0 + T_amb0*0.5*sin(pi*i*timestepLength/24);

    %solve the system of nonlinear equations
    %with initial guess of Tw, Tg, Tb
    fun = @(x) SolarEnergyBalances(x, G, Tamb, wLevel,A); 
    x0 = [20,21,22];    
    x = fsolve(fun,x0); 

    %from solution, use Tw to find evaporative heat transfer, and
    %from that, hourly water output
    Pw = exp(25.317-5144/(273.15+x(2)));
    Pg = exp(25.317-5144/(273.15+x(1)));
    Hc_gw = 0.884*(x(2)-x(1) + x(2)*(Pw-Pg)/(268.9*1000 - Pw))^(1/3);
    He_gw = (16.273*10^(-3))*Hc_gw*(Pw-Pg)/(x(2)-x(1));
    Litres_hourly_clean_water = A*He_gw*(x(2)-x(1))*(3600*timestepLength)/latentHeatVaporization;

    Daily_Output2([i,1])=Litres_hourly_clean_water;
    Efficiency_Matrix([i,1]) = A*He_gw*(x(2)-x(1))/G;
end

%Daily_Output = zeros([amountSteps,1]);
Daily_Output2([1,1])=0;

total_G_KWH = sum(Irradiation_Matrix)
total_daily_water_output1 = sum(Daily_Output2)


for i = 1:amountSteps
    %Savethis_Matrix([i,1]) = Daily_Output([i,1]);
    %Savethis_Matrix([i,1]);
end

%fprintf('\n\n Solar still daily water output = %1.2f L \n\n',total_daily_water_output1);
plot(Time_Matrix,Daily_Output2/timestepLength,'r--o','LineWidth',0.7)
%title('Potable Water Output Over Time for Single Basin Solar Still')
%xlabel('Time since sunrise [hrs]')
%ylabel('Potable water output [L/hr]')
%xlim([0 12.2])
%ylim([-0.1 1])
%legend('Best Net Output = 6.29L','Worst Net Output = 4.13L')

%%

%clear all
%close all
%clc

timestepLength = 0.25; %steplength of one is 1 hour long, 0.5 is 30 min etc..
amountSteps = 24/timestepLength;    % total amount of time steps required 

Daily_Output3 = zeros([amountSteps,1]);
Time_Matrix = zeros([amountSteps,1]);
Efficiency_Matrix = zeros([amountSteps,1]);
Irradiation_Matrix = zeros([amountSteps,1]);
%Savethis_Matrix = zeros([amountSteps,1]);


for i = 1:amountSteps
    Time_Matrix([i,1])= i*timestepLength-1;
    Time_Matrix([1,1])=0;
end


%initialize variables that cannot be included in energy balance
Gmax = 460;  %max solar irradiance in worst conditions
%Gmax = 890;  %max solar irradiance in worst conditions
Gmax = (890+460)/2;

latentHeatVaporization = 2.25*10^6;
T_amb0 = 20;                    %initial ambient temperature
startingWaterLevel = 0.02;      %initial water level
A = 1.5;                  %water and basin surface area (only for water level changes in this code)


for i = 1:amountSteps
    if i == 1
        wLevel = startingWaterLevel;
    else
        wLevel = wLevel - (Litres_hourly_clean_water/A)/1000;
    end
    
    %Solar irradiation as a function of time (hrs)
    G = (Gmax/2)*(sin(pi*(i*timestepLength-1)/9.83))+abs((Gmax/2)*(sin(pi*(i*timestepLength-1)/9.83)));
    Irradiation_Matrix([i,1]) = G*timestepLength;
    
    Tamb = T_amb0 + T_amb0*0.5*sin(pi*i*timestepLength/24);

    %solve the system of nonlinear equations
    %with initial guess of Tw, Tg, Tb
    fun = @(x) SolarEnergyBalances(x, G, Tamb, wLevel,A); 
    x0 = [20,21,22];    
    x = fsolve(fun,x0); 

    %from solution, use Tw to find evaporative heat transfer, and
    %from that, hourly water output
    Pw = exp(25.317-5144/(273.15+x(2)));
    Pg = exp(25.317-5144/(273.15+x(1)));
    Hc_gw = 0.884*(x(2)-x(1) + x(2)*(Pw-Pg)/(268.9*1000 - Pw))^(1/3);
    He_gw = (16.273*10^(-3))*Hc_gw*(Pw-Pg)/(x(2)-x(1));
    Litres_hourly_clean_water = A*He_gw*(x(2)-x(1))*(3600*timestepLength)/latentHeatVaporization;

    Daily_Output3([i,1])=Litres_hourly_clean_water;
    Efficiency_Matrix([i,1]) = A*He_gw*(x(2)-x(1))/G;
end

%Daily_Output = zeros([amountSteps,1]);
Daily_Output3([1,1])=0;

total_G_KWH = sum(Irradiation_Matrix)
total_daily_water_output1 = sum(Daily_Output3)


for i = 1:amountSteps
    %Savethis_Matrix([i,1]) = Daily_Output([i,1]);
    %Savethis_Matrix([i,1]);
end

%fprintf('\n\n Solar still daily water output = %1.2f L \n\n',total_daily_water_output1);
%plot(Time_Matrix,Daily_Output3/timestepLength,'r--o','LineWidth',0.7)
%title('Potable Water Output Over Time for Single Basin Solar Still')
%xlabel('Time since sunrise [hrs]')
%ylabel('Potable water output [L/hr]')
%xlim([0 12.2])
%ylim([-0.1 1])
%legend('Best Net Output = 6.29L','Worst Net Output = 4.13L')
%%

figure;
pause(3);
for i = 1:amountSteps
    hold on;
    plot(Time_Matrix([i,1]),Daily_Output([i,1])/timestepLength,'ro','LineWidth',1.5)
    plot(Time_Matrix([i,1]),Daily_Output1([i,1])/timestepLength,'go','LineWidth',1.5)
    plot(Time_Matrix([i,1]),Daily_Output2([i,1])/timestepLength,'bo','LineWidth',1.5)
    plot(Time_Matrix([i,1]),Daily_Output3([i,1])/timestepLength,'k*','LineWidth',1.5)
    legend('Wind velocity 1m/s -> 4m/s = 5.13 L','Insulation thickness 0.02m -> 0.01m = 2.86 L', 'Basin absorbtivity 0.87 -> 0.7 = 3.51 L', 'Reference output = 4.55 L');
    xlim([0 12.2])
    ylim([-0.1 1])
    title('Potable Water Output Over Time for Single Basin Solar Still')
    xlabel('Time since sunrise [hrs]')
    ylabel('Potable water output [L/hr]')
    drawnow;
end