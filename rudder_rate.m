%% rudder design

b = 1;              % wingspan [m]
b_ai = b*0.1;     % inner coordinates of rudder [m]
b_ao = b*0.9;     % outer coordinates of rudder [m]
Ca_C = 0.3;         % ratio of aileron-to-vertical tail wing surface
da_max = 30*pi/180; % maximum rudder deflection [rad]

Ct = 0.16;   % wing tip chord [m]
Cr = 0.2;   % wing root chord [m]
lambda = Ct / Cr;   % tip, root chord ratio (=taper ratio)
Sv = 0.972/2;          % vertical tail wing surface [m^2]
S = Sv;    % aircraft vertical tail wing surface
I_zz = 0.4158;       % inertia of z axis
CSB = 2.1;         % tail wing sideforce effectiveness [/rad]
Ta = -6.624*(Ca_C)^4+12.07*(Ca_C)^3-8.292*(Ca_C)^2+3.295*(Ca_C)+0.004942;   % rudder's efficiency
Clda = 2*Cr*CSB*Ta/Sv/b*(((b_ao-b_ai)^2)/2+2/3*(lambda-1)/b*((b_ao-b_ai)^3));   % rudder yawing effectiveness [/rad]
CYr = -0.24691;     % yaw damping [/rad]
Cla = Clda*da_max;  % rudder yaw effectiveness
V = 17; % velocity [m/s]
p = 1.225;  % density of air [kg/m^3]
q = p/2*V^2; % [kg/m/s^2]
L_A = q*Sv*Cla*b;   % moment of yaw at maximum rudder deflection
p_ss = -Clda/CYr*da_max*2*V/b; % steady state roll rate [rad/s]
phi_1 = p_ss^2/2/L_A*I_zz*log(p_ss^2);  % bank angle [rad]
p_dot = p_ss^2/2/phi_1; % the rate of roll [rad/s^2]
phi_des = 30*pi/180;    % [rad]
if phi_1>phi_des
    t2 = sqrt(2*phi_des/p_dot); 
else 
    t2 = sqrt(2*phi_1/p_dot)+(phi_des-phi_1)/p_ss;
end
% t2값이  1.3s 보다 작아야 좋은 설계