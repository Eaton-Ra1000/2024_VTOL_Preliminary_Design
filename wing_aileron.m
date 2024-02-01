clear all
close all
clc

%% main wing airfoil is NACA6412

%% 3-D wing condition and definition

b = 2;   % wing span [m]
Ct = 0.2;   % wing tip chord [m]
Cr = 0.26;   % wing root chord [m]
c = @(y) Cr+(Ct-Cr)/(b/2)*y;    % function for chord at y
Sw = 2*integral(c,0,b/2);  % wing area [m^2]
lambda = Ct / Cr;   % tip, root chord ratio (=taper ratio)
A = (b ^ 2) / Sw;   % aspect ratio
ALE = 0;    % leading edge sweep angle [rad]
epsilonr = 0;   % root 뒷틀림각
epsilont = 0;   % tip 뒷틀림각
diangle = 0;    % dihedral angle [rad]

%% center point of chord sweep angle: 모든 시위 중앙을 지나는 선의 뒷젖침각

centerC = Cr/2 - tan(ALE)*b/2 - (Ct / 2);    % line along center point of chord [m]
Ac2 = atan(centerC / (b / 2));  % center point of chord sweep angle [rad]

%% wing lift effectiveness: 양력선 기울기

Minf = 0.058;   % mach 0.05, 17m/s 
Beta = sqrt(1 - (Minf ^ 2));    % frandtl subsonic compressivity coefficient    양력선 기울기
cla = 0.105*180/pi;  % section lift effectiveness [/rad]
kappa = cla / (2 * pi); % airfoil camber degree dimentionless coefficient [/rad]
CLa = (2 * pi * A) / (2 + sqrt(((((A ^ 2) * (Beta ^ 2)) / (kappa ^ 2)) * (1 + ((tan(Ac2) ^ 2) / Beta ^ 2)) + 4)));    % wing lift effectiveness [/rad]  주익 양력 효과

%% wing zero-lift angle of attack: 날개에서 양력이 발생하지 않는 특정 받음각

a_0 =  -6.18*pi/180;    % 2D airfoil의 양력이 발생하지 않는 받음각 [rad]     %%%% CFD로 구해야함
epsilon = @(y) epsilonr+(epsilont-epsilonr)/(b/2)*y;    % y에서 뒤틀림 각 [rad]
fun_a_0wing = @(y) (a_0 - epsilon(y)).*c(y);    
a_0wing = 2/Sw*integral(fun_a_0wing,0,b/2);     % 3d에서 양력이 발생하지 않는 받음각 [rad]

%% 평균 공력 시위

x_LE = @(y) tan(ALE).*y;    % LE x 좌표 [m]
fun_c_bar = @(y) c(y).^2;
c_bar = 2/Sw*integral(fun_c_bar,0,b/2);     % 평균 공력시위 [m]
fun_X_LE_MAC = @(y) x_LE(y).*c(y);
X_LE_MAC = 2/Sw*integral(fun_X_LE_MAC,0,b/2);   % LE ~ MAC 까지 거리
X_AC_wing = X_LE_MAC + 0.25*c_bar ;    % 공력중심 x좌표 [m]
fun_Y_MAC = @(y) y.*c(y);
Y_MAC = 2/Sw*integral(fun_Y_MAC,0,b/2); % 공력중심 y좌표 [m]

%% 날개 피칭 모멘트 계수

c_m_ac = -0.14;     %%%%%% 공력 중심에서의 피칭 모멘트 (양력이 0일때 모멘트 계수 c_m)  %%%% CFD 활용
fun_C_M_AC1 = @(y) c_m_ac.*c(y).^2;
x_ac = @(y) x_LE(y)+0.25.*c(y);
fun_C_M_AC2 = @(y) cla.*(a_0wing+epsilon(y)-a_0).*(x_ac(y)-X_AC_wing).*c(y);
C_M_AC = 2/Sw/c_bar*(integral(fun_C_M_AC1,0,b/2)-integral(fun_C_M_AC2,0,b/2));  % wing pitching moment codfficient 날개 피칭 모멘트 계수

%% 롤링 모멘트 계수 (받음각이 a_wing, side slip이 B일 때)

B = 2*pi/180;  % side slip 옆 미끄러짐각 [rad]
a_wing = 2*pi/180;     % 받음각 [rad]

C_L_roll_ALE = CLa*(a_wing - a_0wing)*Y_MAC/b*(-sin(2*ALE)*sin(2*B));    % 뒷젖침각과 옆 미끄러짐각으로 생기는 롤링 모멘트 계수
C_L_roll_diangle = -diangle*CLa*B*Y_MAC/b;   % 상반각과 옆 미끄러워짐각으로 생기는 롤링 모멘트 계수
C_L_roll = C_L_roll_ALE+C_L_roll_diangle;    % 주익의 롤링 모멘트 계수 [단위없음]

%% 항력계수 (받음각 a_wing일 때)

t_c = 0.12;  % thickness ratio of the airfoil section
S_wet = 2*Sw;  % 공기에 노출된 날개 면적 (2*Sw와 비슷)      %%% 모델링 이후 계산가능
p = 1.225;  % density of air [kg/m^3]
mu = 1.54*10^(-5); % 공기 점성계수 [m^2/s]
V = 17;     % 속도    [m/s]
l = c(Y_MAC);   % 특성길이
R_l = p*V*l/mu; % Reynolds number

C_f = 0.0053632; % skin friction coefficient 표면 마찰계수    (Figure 5.31) (aerotoolbox사이트 계산기 있음)
C_D_P = C_f*(1+2*(t_c)+100*(t_c)^4)*S_wet/Sw;   % 마찰 항력계수
C_L = CLa*(a_wing-a_0wing);     % 양력계수
e = 0.86;       % osweld's efficiency factor    (Figure osweld's effiency factor)
C_D_I = C_L^2/pi/A/e;     % 유도 항력계수
C_D = C_D_P+C_D_I;   % 항력계수 [/rad]

%% Aileron rolling effectiveness : 에일러론 롤링 효과

b_ai = b/2*0.5;  % inner coordinates of aileron [m]
b_ao = b/2*0.7;  % outer coordinates of aileron [m]
Ca_C = 0.2;  % ratio of aileron chord to wing chord
da_max = 20*pi/180; % maximum aileron deflection [rad]

cl_delta_theory = 3.7; % Theoretical section-lift effectiveness [/rad]  (Figure 5.14)
cla_theory = 2*pi;     % theory lift effecriveness [/rad]
cla_over_cla_theory = cla/cla_theory;   % ratio of theory 
cl_over_cl_theory = 0.91; % Empirical section-lift effectiveness  (Figure 5.15)

C_l_delta_aileron = 1/Beta*(cl_over_cl_theory)*cl_delta_theory;
fun_C_L_delta_aileron = @(y) C_l_delta_aileron.*c(y).*y;
C_L_delta_aileron = 2/Sw/b*integral(fun_C_L_delta_aileron,b_ai,b_ao);   % 에일러론으로 롤링 모멘트 효과 [/rad]

%% Aileron yawing moment effectiveness : 에일러론 요잉 모멘트 효과

del_delta = 20*pi/180;     % 에일러론 변화 각 [rad]
del_c_d = 0.02;     % delta 각에 따른 항력 계수 (Figure 5.17)

fun_C_N_delta_aileron = @(y) del_c_d/del_delta.*c(y).*y;
C_N_delta_aileron = -2/Sw/b*integral(fun_C_N_delta_aileron,b_ai,b_ao);  % 에일러론 요잉 모멘트 효과 [/rad]
