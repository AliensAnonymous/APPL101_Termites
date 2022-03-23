clear; clc; close all;

% Starting Values
% Simple starting model so will not include thermal diffusivity on thermal
% conductivity
% This is not accurate numericaly, more just a proof of concept with
% extremely simplified equations
% https://www.researchgate.net/publication/304719050_An_Investigation_into_the_Thermal_Properties_of_Termite_Mound_Clay_Applicable_to_Grain_Silo_Construction#pf3

% Specific Heat of Air - J/kg/C
airSH = 1003.5;

% Specific Heat of TMC (Termite Mound Clay) - J/kg/C 
tmcSH = 2576.9;

% Thermal conductivity of TMC - W/mK
tmcTC = 0.21;

% Thermal conductivity of Air - W/mK
airTC = 0.025;

% Density of Air - kg/m3
airDensity = 1.225;

% Density of TMC - kg/m3
tmcDensity = 1833;

% Air Material Resistance
airResist = 1/(airSH * airDensity);

% TMC Material Resistance
tmcResist = 1/(tmcSH * tmcDensity);

simWidth = 40;
simHeight = 40;
simLast = zeros(simWidth, simHeight);
sim = zeros(simWidth, simHeight);
simTimeSteps = 1000;
heatSource = 38;
coldSource = 14;
totalTimeLength = 10000000;
resistMap = zeros(simWidth, simHeight);
    % Material map values will be thermal resistance of respective
    % materials
    
% Air Adjusted Conductivity
airCon = airTC * (totalTimeLength/simTimeSteps);

% TMC Adjusted Conductivity
tmcCon = tmcTC * (totalTimeLength/simTimeSteps);
 
conductMap = zeros(simWidth, simHeight);

f1 = figure('Position',[1000, 100, 500, 500]);

for y = 1:simHeight
    if y <= simHeight/2
        resistMap(y, :) = airResist;
        conductMap(y,:) = airCon;
    else
        resistMap(y,:) = tmcResist;
        conductMap(y,:) = tmcCon;
    end 
    
    if y == 1
        simLast(1,:) = heatSource;
    elseif y == simHeight
        simLast(simHeight,:) = coldSource;
    else
        simLast(y,:) = (heatSource + coldSource)/2;
    end
end
sim = simLast;

% Start Calcuations

for t = 1:simTimeSteps
    for i = 2:simHeight-1
        for j= 1:simWidth
            if j == 1
                sim(i,j) = simLast(i,j) + resistMap(i,j) * ((conductMap(i - 1,j) * (simLast(i - 1,j) - simLast(i,j))) + (conductMap(i + 1,j) * (simLast(i + 1,j) - simLast(i,j))) + (conductMap(i,j + 1) * (simLast(i,j + 1) - simLast(i,j))));
            elseif j == simWidth
                sim(i,j) = simLast(i,j) + resistMap(i,j) * ((conductMap(i - 1,j) * (simLast(i - 1,j) - simLast(i,j))) + (conductMap(i + 1,j) * (simLast(i + 1,j) - simLast(i,j))) + (conductMap(i,j - 1) * (simLast(i,j - 1) - simLast(i,j))));
            else 
                sim(i,j) = simLast(i,j) + resistMap(i,j) * ((conductMap(i - 1,j) * (simLast(i - 1,j) - simLast(i,j))) + (conductMap(i + 1,j) * (simLast(i + 1,j) - simLast(i,j))) + (conductMap(i,j - 1) * (simLast(i,j - 1) - simLast(i,j))) + (conductMap(i,j + 1) * (simLast(i,j + 1) - simLast(i,j))));
            end
        end
    end  
    simLast = sim;
    imagesc(sim, [14 38]);
    drawnow;
end
    