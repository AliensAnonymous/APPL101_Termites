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

% Cell size - m
size = 1;

% Air Material Resistance
airResist = (1/(airSH * airDensity))*(size^2);

% TMC Material Resistance
tmcResist = (1/(tmcSH * tmcDensity))*(size&2);

% Incoming insolation - W/m2
insolation = 1380;


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

% Setting starting values, material positions
for y = 1:simHeight
    if y <= simHeight/2
        resistMap(y, :) = airResist;
        conductMap(y,:) = airCon;
        
        % Uncomment to change starting conditions
        %simLast(y,:) = heatSource;
    else
        resistMap(y,:) = tmcResist;
        conductMap(y,:) = tmcCon;
        
        % Uncomment to change starting conditions
        %simLast(y,:) = coldSource;
    end 
    
    % Setting boundary conditions
     if y == 1
         simLast(1,:) = heatSource;
     elseif y == simHeight
         simLast(simHeight,:) = coldSource;
     % Setting background starting temperature as halfway between hot and
     % cold
     else
         simLast(y,:) = (heatSource + coldSource)/2;
     end
end
sim = simLast;

% Start Calcuations
top = 0;

for t = 1:simTimeSteps
    for j= 1:simWidth
        for i = 2:simHeight-1
            % Boundary conditions
            if j == 1
                sim(i,j) = simLast(i,j) + resistMap(i,j) * ((conductMap(i - 1,j) * (simLast(i - 1,j) - simLast(i,j))) + (conductMap(i + 1,j) * (simLast(i + 1,j) - simLast(i,j))) + (conductMap(i,j + 1) * (simLast(i,j + 1) - simLast(i,j))));
            elseif j == simWidth
                sim(i,j) = simLast(i,j) + resistMap(i,j) * ((conductMap(i - 1,j) * (simLast(i - 1,j) - simLast(i,j))) + (conductMap(i + 1,j) * (simLast(i + 1,j) - simLast(i,j))) + (conductMap(i,j - 1) * (simLast(i,j - 1) - simLast(i,j))));
            % Normal conditions
            else 
                sim(i,j) = simLast(i,j) + resistMap(i,j) * ((conductMap(i - 1,j) * (simLast(i - 1,j) - simLast(i,j))) + (conductMap(i + 1,j) * (simLast(i + 1,j) - simLast(i,j))) + (conductMap(i,j - 1) * (simLast(i,j - 1) - simLast(i,j))) + (conductMap(i,j + 1) * (simLast(i,j + 1) - simLast(i,j))));
            end
            % Adding insolation on the upper layer of TMC 
            if top == 0 && resistMap(i,j) == tmcResist
                top = 1;
                sim(i,j) = sim(i,j) + resistMap(i,j) * size * (insolation * totalTimeLength/simTimeSteps*size);
            end
        end
        top = 0;
    end  
    % Updating and drawing sim
    simLast = sim;
    imagesc(sim, [14 45]);
    drawnow;
end
    