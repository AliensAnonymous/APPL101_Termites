clear; clc; close all;
totalTimeLength = 1000;  
simTimeSteps = 1000;
% https://www.researchgate.net/publication/304719050_An_Investigation_into_the_Thermal_Properties_of_Termite_Mound_Clay_Applicable_to_Grain_Silo_Construction#pf3
airSH = 1003.5;% Specific Heat of Air - J/kg/C
tmcSH = 2576.9;% Specific Heat of TMC (Termite Mound Clay) - J/kg/C 
tmcTC = 0.21;% Thermal conductivity of TMC - W/mK
airTC = 0.025;% Thermal conductivity of Air - W/mK
airDensity = 1.225;% Density of Air - kg/m3
tmcDensity = 1833;% Density of TMC - kg/m3
airResist = 1/(airSH * airDensity);% Air Material Resistance
tmcResist = 1/(tmcSH * tmcDensity);% TMC Material Resistance
airCon = airTC * (totalTimeLength/simTimeSteps); % Air Adjusted Conductivity
tmcCon = tmcTC * (totalTimeLength/simTimeSteps); % TMC Adjusted Conductivity
h = 5;
b = 100;
simX = 301;
simY = 475;
temp = zeros(simY, simX);
temp1 = zeros(simY, simX);
heatSource = 38;
coldSource = 14;


image = imread('moundImage.jpg');
image = imcrop(image,[75,0,300,475]);

map = image ~= 255;

for i = 2:(simY-1)
    for j = 1:simX
        if map(i,j) == 1
            temp(i,j) = coldSource;
        else
            temp(i,j) = (heatSource + coldSource)/2;
        end
    end
end

%map = ones(simY, simX);
%simXhalf = int64(simX/2)
%simYhalf = int64(simY/2)
%map(1:simXhalf, :) = 0;
%map(simYhalf:simY, :) = 1;
%temp(2:simYhalf, :) = (heatSource + coldSource)/2;
%temp(simYhalf+1:simY, :) = coldSource;
temp(1, :) = heatSource;
temp(simY, :) = coldSource;
dt = totalTimeLength/simTimeSteps;
dx = .01;

f1 = figure('Position',[1000, 100, 602, 950]);
temp1 = temp;
imagesc(temp1, [14 38]);
drawnow;

dTemp = 0;
q_right = 0;
q_down = 0;
q_up = zeros(simX, 1);
% Start Calcuations
for t = 1:simTimeSteps
    %calculate heat for first row
    for x = 1:simX
        q_up(x) = conduct(airTC, temp(1, x), temp(2, x), dt, dx);
    end
    %calculate heat for all other rows
    for y = 2:simY - 1
        q_left = 0;
        q_up = -1*q_up;
        for x = 1:simX-1
           if map(y, x) == 0
               k = airTC;
               c = airSH;
               p = airDensity;
           else
               k = tmcTC;
               c = tmcSH;
               p = tmcDensity;
           end
           q0 = q_up(x) + q_left;
           [dTemp, q_right, q_down] = heat(q0, dt, dx, temp(y, x), temp(y, x+1), temp(y+1, x), k, p, c, map(y, x), map(y, x+1), map(y+1, x));
           q_left = -q_right;
           q_up(x) = q_down;
           temp1(y, x) = temp(y, x) + dTemp;
        end
        if (map(y, simX) == map(y+1, simX))
            q_down = conduct(k, temp(y, simX), temp(y+1, simX), dt, dx);
        else
            q_down = convect(k, temp(y, simX), temp(y+1, simX), dt, dx);
        end
        q_total = q_left + q_down + q_up(simX);
        q_up(simX) = q_down;
        m = dx^3*p;
        dTemp = q2T(q_total, m, c);
        temp1(y, simX) = temp(y, simX) + dTemp;
    end  
    temp = temp1;
    imagesc(temp1, [14 38]);
    drawnow;
end

function [dT, q_right, q_down] = heat(q0, dt, dx, T1, T2, T3, k, p, c, m_self, m_right, m_down)
    %calculate right cell
    if (m_self == m_right)
        q_right = conduct(k, T1, T2, dt, dx);
    else
        q_right = convect(k, T1, T2, dt, dx);
    end
    %calculate bottom cell
    if (m_self == m_down)
        q_down = conduct(k, T1, T3, dt, dx);
    else
        q_down = convect(k, T1, T3, dt, dx);
    end
    q_total = q_down + q_right + q0;
    m = dx^3*p;
    dT = q2T(q_total, m, c);
end

function dTemp = q2T(q, m, c)
    dTemp = -q/(m*c);
end

function q = conduct(k, T1, T2, dt, dx)
    q = k*dx*(T1-T2)*dt;
end
    
function q = convect(k, T1, T2, dt, dx)
    q = k*dx^2*(T1-T2)*dt;
end