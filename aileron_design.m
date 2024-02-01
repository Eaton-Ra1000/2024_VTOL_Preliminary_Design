%% aileron design

clear all
close all
clc

b = 2;              % wingspan [m]
b_ai = b/2*0.55;     % inner coordinates of aileron [m]
b_ao = b/2*0.9;     % outer coordinates of aileron [m]
Sa_S = 0.0875;         % ratio of aileron-to-wing surface
da_max = 25*pi/180; % maximum aileron deflection [rad]

Ct = 0.2;   % wing tip chord [m]
Cr = 0.26;   % wing root chord [m]
lambda = Ct / Cr;   % tip, root chord ratio (=taper ratio)
Sw = 0.46;          % wing surface [m^2]
S = 0.46+0.0972;    % aircraft wing surface
I_xx = 0.375;       % inertia of x axis
CLa = 4.84;         % wing lift effectiveness [/rad]
Ta = -6.624*(Sa_S)^4+12.07*(Sa_S)^3-8.292*(Sa_S)^2+3.295*(Sa_S)+0.004942;   % aileron's efficiency
Clda = 2*Cr*CLa*Ta/Sw/b*(((b_ao-b_ai)^2)/2+2/3*(lambda-1)/b*((b_ao-b_ai)^3));   % aileron roll effectiveness [/rad]
Clp = -0.54966;     % roll damping [/rad]
Cla = Clda*da_max;  % aileron roll effectiveness
V = 17; % velocity [m/s]
p = 1.225;  % density of air [kg/m^3]
q = p/2*V^2; % [kg/m/s^2]
L_A = q*Sw*Cla*b;   % moment of roll at maximum aileron deflection
p_ss = -Clda/Clp*da_max*2*V/b; % steady state roll rate [rad/s]
phi_1 = p_ss^2/2/L_A*I_xx*log(p_ss^2);  % bank angle [rad]
p_dot = p_ss^2/2/phi_1; % the rate of roll [rad/s^2]
phi_des = 40*pi/180;    % [rad]
if phi_1>phi_des
    t2 = sqrt(2*phi_des/p_dot); 
else 
    t2 = sqrt(2*phi_1/p_dot)+(phi_des-phi_1)/p_ss;
end
% t2값이  1.7s 보다 작아야 좋은 설계