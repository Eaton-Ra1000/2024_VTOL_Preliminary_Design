clear all
close all
clc
%% main wing airfoil is NACA6412


%% 3-D wing condition and definition

b = 2;   % wing span [m]
Sw = 0.46;  % wing area [m^2]
Ct = 0.2;   % wing tip chord [m]
Cr = 0.26;   % wing root chord [m]
lambda = Ct / Cr;   % tip, root chord ratio (=taper ratio)
A = (b ^ 2) / Sw;   % aspect ratio
ALE = 0;    % leading edge sweep angle [degree]
c(y) = Cr+(Ct-Cr)/(b/2)*y;  % function for chord at y
epsilonr = 2;   %%%%%% root 뒷틀림각, 미정
epsilont = -1;   %%%%% tip 뒷틀림각, 미정
diangle = 7.5;    %%%%% dihedral angle [degree]

%% center point of chord sweep angle: 모든 시위 중앙을 지나는 선의 뒷젖침각

centerC = Cr/2 - (Ct / 2);    % line along center point of chord [m]
Ac2 = atan(centerC / (b / 2));  % center point of chord sweep angle [rad]

%% wing lift effectiveness: 양력선 기울기

Minf = 0.058;   % mach 0.05, 17m/s 
Beta = sqrt(1 - (Minf ^ 2));    % frandtl subsonic compressivity coefficient
cla = 0.105*180/pi;  % section lift effectiveness [/rad]
kappa = cla / (2 * pi); % airfoil camber degree dimentionless coefficient [/degree]

CLa = (2 * pi * A) / (2 + sqrt(((((A ^ 2) * (Beta ^ 2)) / (kappa ^ 2)) * (1 + ((tan(Ac2) ^ 2) / Beta ^ 2)) + 4)));    % wing lift effectiveness [/rad]
%% wing zero-lift angle of attack: 날개에서 양력이 발생하지 않는 특정 받음각

a_0 =  -6.2;    %%%%%% 2D airfoil의 양력이 발생하지 않는 받음각_NACA6412@RE=200,000
e(y) = epsilonr+(epsilont-epsilonr)/(b/2)*y;    % y에서 뒤틀림 각
a_0wing = 2/Sw*int((a_0-e(y))*c(y),0,b/2);      % 3차원 날개에서 양력이 발생하지 않는 받음각
c_bar = 2/Sw*int(c(y)^2,0,b/2);                 %평균공력시위 길이
x_LE(y) = 0;
X_LE_MAC = 2/Sw*int(x_LE(y)*c(y),0,b/2);        % Y축에서부터 X축길이로 평균 공력 시위 (MAC)가 위치한 곳의 앞전(Leading edge)까지의 거리
Y_MAC = 2/184*int(y*c(y),0,b/2);
X_AC_wing = X_LE_MAC + 0.25*c_bar;
%% wing pitching moment coefficient: 날개 피칭 모멘트 계수

C_M_AC = 2/Sw/c_bar*(int(c_m_ac(y)*c(y)^2,0,b/2)-int(cla(y)*(a_0wing+epsilon(y)-a_0(y)*()),0,b/2)-










%% aileron

b_ai = 0.5;  % inner coordinates of aileron [m]
b_ao = 0.7;  % outer coordinates of aileron [m]
Ca_C = 0.2;  % ratio of aileron chord to wing chord
da_max = 20*pi/180; % maximum aileron deflection [rad]
Ta = -6.624*(Ca_C)^4+12.07*(Ca_C)^3-8.292*(Ca_C)^2+3.295*(Ca_C)+0.004942;
Clda = 2*Cr*CLa*Ta/Sw/b*((b_ao-b_ai)^2/2+2/3*(lambda-1)/b*(b_ao-b_ai)^3);   % lateral control derivative of the aileron
Cla = Clda*da_max;  % aerodynamic coefficient of the aileron
Vs = 17;    %%%%%%%%%
Vt = Vs*1.3;    %%%%%%
p = 1.225;  % density of air [kg/m^3]
q = p/2*Vt^2; %%%%%%%%
L_A = q*Sw*Cla*b;   % moment of roll at maximum aileron deflection
C_DR = 0.9;
Y_d = 0.4;  % coordinate of mean drag of three lifting surfaces(m)
Sh = 0.1;   % horizon tail wing area [m^2]
Sv = 0.02;  % vertical tail wing area [m^2]
P_ss = sqrt(2*L_A/p/(Sw+Sh+Sv)/C_DR/Y_d^3) % Steady State Roll Rate [rad/s]



