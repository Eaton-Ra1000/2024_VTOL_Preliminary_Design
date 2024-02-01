% elevator design

clear all
close all
clc

de = -20*pi/180; % elevator deflection [rad]
V_a = 17;   % velecity [m/s]
c = 0.2313; % charateristic length [m]
Sw = 0.46;  % wing surface [m^2]
Sh = 0.0972;    % horizon tail surface [m^2]
C_M0F = -0.0288; %%%% fuselage moment coeff 
C_M_AC = -0.14; % wing moment coeff
C_L_a = 4.84;   % wing lift coeff
C_L_ah = 3.38;  % horizon tail lift coeff
C_m_q = -12.65643;  % pitch damping
p = 1.225;  % density of air [kg/m^3]
a_0wing = -0.10786; % zero lift aoa [rad]
a_wing = 0*pi/180;  % wing aoa [rad]
downwash_a_w = 0.43;    % 수평 꼬리날개 다운워시 각 계수
e0 = 0.0382;    % downwash when aoa=0
downwash = e0+downwash_a_w*(a_wing-a_0wing);   % 수평 꼬리날개 다운워시 각 [rad]
a_w = a_wing+downwash;   % wing aoa [rad]
a_h = a_wing-downwash;  % tail wing aoa [rad]
I_yy = 0.106134;    % inertia y axis [kg*m^2]
q_dot = 15*pi/180;     % 설정 값 [rad/s^2]
M_0F = C_M0F*p*Sw*c*V_a^2/2; % fuselage moment
M_ac_w = C_M_AC*p*Sw*c*V_a^2/2; % wing moment
L_w = C_L_a*a_w*(0.003)*0.9*0.211;  % wing lift
M_q = 1/4*p*V_a*Sw*c^2*C_m_q;   % moment of pitch velocity

L_h = (I_yy*q_dot-M_0F-L_w*(0.003)-M_ac_w-M_q)/(-0.679);    % tail wing lift
C_L_h = 2*L_h/p/V_a^2/Sh;   % tail wing coeff
Te_des = (C_L_h-C_L_ah*a_h)/C_L_ah/de; % desire elevator effectiveness

syms Ca_C f
f = -6.624.*(Ca_C).^4+12.07.*(Ca_C).^3-8.292.*(Ca_C).^2+3.295.*(Ca_C)+0.004942-Te_des;
Ca_C_des = vpasolve(f,Ca_C,[0,1])   % desire chord ratio

C_m0 = -0.03553;
C_L0 = 0.45132;
C_m_a = 0.014;
C_m_da = -0.5125;
C_L_da = 3.598;
de_cruse = (C_m0*C_L_a+(-C_L0)*C_m_a)/(C_L_a*C_m_da-C_m_a*C_L_da)
a_h_stall = 5.3*pi/180;
a_to = a_h_stall+downwash

