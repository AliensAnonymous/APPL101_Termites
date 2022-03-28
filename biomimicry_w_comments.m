clear; clc; close all;
%Initializing Initial Conditions
totalTimeLength = 1000;  
simTimeSteps = 1000;
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
simX = 301; %x dimension of sim grid (adjusted to termite mound image)
simY = 475; %y dimension of sim grid (adjusted to termite mound image)

%initialize temperature arrays
temp = zeros(simY, simX);
temp1 = zeros(simY, simX);

heatSource = 38; %standard max. air temperature 
coldSource = 14; %standard ground temperature

%import termite mound image
image = imread('moundImage.jpeg');
image = imcrop(image,[75,0,300,475]);

%overlay map array onto mound image
map = image ~= 255;

%initialize map array with initial temperature conditions for corresponding
%material types
for i = 2:(simY-1)
    for j = 1:simX
        if map(i,j) == 1
            temp(i,j) = coldSource; 
        else
            temp(i,j) = (heatSource + coldSource)/2;
        end
    end
end

temp(1, :) = heatSource; %top of temperature array is max heat (source)
temp(simY, :) = coldSource; %tmc/ground half of map is uniformly coldest (=cold source)

dt = totalTimeLength/simTimeSteps; %change in time over simulation
dx = .01;

%draws initial heat map in temperature gradient
f1 = figure('Position',[1000, 100, 602, 950]);
temp1 = temp;
imagesc(temp1, [14 38]);
drawnow;

dTemp = 0; %initialize change in temperature

%initialize heat of neighbors
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
           if map(y, x) == 0 %setting initial conditions for air 
               k = airTC;
               c = airSH;
               p = airDensity;
           else %setting initial conditions for tmc
               k = tmcTC;
               c = tmcSH;
               p = tmcDensity;
           end
           q0 = q_up(x) + q_left; 
           [dTemp, q_right, q_down] = heat(q0, dt, dx, temp(y, x), temp(y, x+1), temp(y+1, x), k, p, c, map(y, x), map(y, x+1), map(y+1, x));
           q_left = -q_right; %heat of left neighbor is same as right neighbor
           q_up(x) = q_down; %bottom neighbor becomes top neighbor for next iteration
           temp1(y, x) = temp(y, x) + dTemp;
        end
        if (map(y, simX) == map(y+1, simX)) 
            q_down = conduct(k, temp(y, simX), temp(y+1, simX), dt, dx); %account for thermal conductivity
        else
            q_down = convect(k, temp(y, simX), temp(y+1, simX), dt, dx); %account for convective effects of the air inside the mound
        end
        q_total = q_left + q_down + q_up(simX); %calculate total heat of neighbors
        q_up(simX) = q_down;
        m = dx^3*p; 
        dTemp = q2T(q_total, m, c); %calculate change in temp. for next cell
        temp1(y, simX) = temp(y, simX) + dTemp;
    end  
    temp = temp1;
    imagesc(temp1, [14 38]); %continuosly update temperature image
    drawnow;
end

%heat function - calculate the temperature of the cells using conductions
%and convection rates
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

%claculate cell's change in temperature 
function dTemp = q2T(q, m, c)
    dTemp = -q/(m*c);
end

%calculate conduction of cell using heat transfer differential equation
function q = conduct(k, T1, T2, dt, dx)
    q = k*dx*(T1-T2)*dt;
end

%calculate convection of cell using using convection diffusion equation
function q = convect(k, T1, T2, dt, dx)
    q = k*dx^2*(T1-T2)*dt;
end 