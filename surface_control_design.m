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
ALE = 0;    % leading edge sweep angle [rad]
c = @(y) Cr+(Ct-Cr)/(b/2)*y;    % function for chord at y
epsilonr = 0;   % root 뒷틀림각
epsilont = 0;   % tip 뒷틀림각
diangle = 0;    % dihedral angle [rad]

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

a_0 =  -6.95;    % 2D airfoil의 양력이 발생하지 않는 받음각 [degree]     %%%%%% CFD로 구해야함
epsilon = @(y) epsilonr+(epsilont-epsilonr)/(b/2)*y;    % y에서 뒤틀림 각
fun_a_0wing = @(y) (a_0 - epsilon(y)).*c(y);    % a_0wing 계산을 위해 지정한 함수
a_0wing = 2/Sw*integral(fun_a_0wing,0,b/2);     % 3d에서 양력이 발생하지 않는 받음각 [degree]

%% 평균 공력 시위 (*외부 식)

fun_c_bar = @(y) c(y).^2;
c_bar = 2/Sw*integral(fun_c_bar,0,b/2);     % 평균 공력시위
X_AC_wing = 0.25*c_bar ;    % 공력중심 x좌표
fun_Y_MAC = @(y) y.*c(y);
Y_MAC = 2/Sw*integral(fun_Y_MAC,0,b/2); % 공력중심 y좌표

%% 날개 피칭 모멘트 계수

c_m_ac = -0.12;     %%%%%% 공력 중심에서의 피칭 모멘트 (양력이 0일때 모멘트 계수??)
fun_C_M_AC1 = @(y) c_m_ac.*c(y).^2;
fun_C_M_AC2 = @(y) cla.*(a_0wing+epsilon(y)-a_0).*(0.25*c(y)-X_AC_wing).*c(y);
C_M_AC = 2/Sw/c_bar*(integral(fun_C_M_AC1,0,b/2)-integral(fun_C_M_AC2,0,b/2));  % wing pitching moment codfficient 날개 피칭 모멘트 계수

%% 롤링 모멘트 계수 (받음각 a_wing일 때) (* 받음각, 옆 미끄러짐각 따른 Variation 가능)

a_wing = 2;     % 받음각 [degree]
C_L_roll_ALE = CLa*(a_wing - a_0wing)*Y_MAC/2*(-sin(2*ALE)*sin(2*Beta));    % 뒷젖침각과 옆 미끄러짐각으로 생기는 롤링 모멘트 계수 (* Sweep 이나 ALE로 구성되어 있음)
C_L_roll_diangle = -diangle*CLa*Beta*Y_MAC/b;   %상반각과 옆 미끄러워짐각으로 생기는 롤링 모멘트 계수
C_L_roll = C_L_roll_ALE+C_L_roll_diangle;

%% 항력계수 (받음각 a_wing일 때) 

t_c = 0.12;  % thickness ratio of the airfoil section
S_wet = 2*Sw;  % 공기에 노출된 날개 면적 (2*Sw와 비슷)
p = 1.225;  % density of air [kg/m^3]
mu = 1.785*10^(-5); % 공기 점성계수
V = 17;     % 속도(* Variant처리)
l = c(Y_MAC);   % 특성길이
R_l = p*V*l/mu; % Reynolds number
C_f = 0.0054629; % skin friction coefficient 표면 마찰계수    (Figure 5.31)
C_D_P = C_f*(1+2*(t_c)+100*(t_c)^4)*S_wet/Sw;   % 마찰 항력계수
C_L = CLa*(a_wing-a_0wing)*pi/180;     % 양력계수
e = 0.84;       % osweld's efficiency factor    (Figure osweld's effiency factor)
C_D_I = C_L^2/pi/A/e;     % 유도 항력계수
C_D = C_D_P+C_D_I;   % 항력계수

%% Aileron lift effectiveness : 에일러론 롤링 효과

b_ai = 0.5;  % inner coordinates of aileron [m]
b_ao = 0.7;  % outer coordinates of aileron [m]
Ca_C = 0.2;  % ratio of aileron chord to wing chord
cl_delta_theory = 3.7; % Theoretical section-lift effectiveness [/rad]  (Figure 5.14)
cla_theory = 2*pi/sqrt(1-Minf^2);   % theory lift effectiveness [/rad] (* 어디서 왔는지 체크)
cl_over_cl_theory = 0.91; % Empirical section-lift effectiveness  (Figure 5.15)

C_l_delta_aileron = 1/Beta*(cl_over_cl_theory)*cl_delta_theory;
fun_C_L_delta_aileron = @(y) C_l_delta_aileron.*c(y).*y;
C_L_delta_aileron = 2/Sw/b*integral(fun_C_L_delta_aileron,b_ai,b_ao);   % 에일러론으로 롤링 모멘트 효과

%% Aileron yawing moment effectiveness : 에일러론 요잉 모멘트 효과

del_delta = 0.349;     % 에일러론 변화 각 [rad]
del_c_d = 0.02;     % delta 각에 따른 항력 계수 (Figure 5.17)

fun_C_N_delta_aileron = @(y) del_c_d/del_delta.*c(y).*y;
C_N_delta_aileron = -2/Sw/b*integral(fun_C_N_delta_aileron,b_ai,b_ao);  % 에일러론 요잉 모멘트 효과





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
P_ss = sqrt(2*L_A/p/(Sw+Sh+Sv)/C_DR/Y_d^3); % Steady State Roll Rate [rad/s]

