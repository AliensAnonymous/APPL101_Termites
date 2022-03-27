clear; clc; close all;
totalTimeLength = 20000;
dt = 0.75;
dx = 0.01;
simTimeSteps = int64(totalTimeLength/dt);

%It should be noted that heats move ONLY by conduction and convection in
%this simulation. In reality, heat rising(that is, air rising as it
%decreases in density) also contributes greatly to the distribution of heat
%in the atmosphere. This air is acting like a solid, and is not moving.


I = imread('termite_mound.png');
[rows, columns, numberOfColorChannels] = size(I);
redChannel = I(:, :, 1);
greenChannel = I(:, :, 2);
blueChannel = I(:, :, 3);
material = (redChannel == 96) & (greenChannel == 56) & (blueChannel == 19);
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
simX = columns;
simY = rows;
temp = zeros(simY, simX);
temp1 = zeros(simY, simX);
heatSource = 38;
coldSource = 14;
map = zeros(simY, simX);
simXhalf = int64(simX/2);
simYhalf = int64(simY/2);
map(material) = 1;
temp(:, :) = (heatSource + coldSource)/2;
temp(material) = 14;
temp(simY, :) = coldSource;
solar = 24*dt*dx^2;

q_x_first = ones(simX);
%set up an initial condition of solar irradiation


f1 = figure('Position',[1000, 100, 500, 500]);
temp1 = temp;
imagesc(temp1, [14 38]);
colorbar;

drawnow;
dTemp = 0;
q_right = 0;
q_down = 0;
q_up = zeros(simX, 1);
t_half = int64(simTimeSteps/2);

% Start Calcuations
%with initial conditions, it acts like the sun just rose and is starting to
%add heat to the system
for t = 1:simTimeSteps
    q_x_first = ones(simX);
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
           if (q_x_first(x) == 1 && map(y+1, x) == 1)
                %a negative heat signifies negative heat-flow, meaning heat is flowing in
                q0 = q0 - solar;  
                q_x_first(x) = 0;
           end
           [dTemp, q_right, q_down] = heat(q0, dt, dx, temp(y, x), temp(y, x+1), temp(y+1, x), k, p, c, map(y, x), map(y, x+1), map(y+1, x));
           q_left = -q_right;
           q_up(x) = q_down;
           temp1(y, x) = temp(y, x) + dTemp;
        end
        %deal with the last row of x (needs to be done seperately)
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
    colorbar
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