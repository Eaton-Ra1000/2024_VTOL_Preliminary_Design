clear all
close all
clc

%% horizon tail wing airfoil is NACA0010

%% 3-D wing condition and definition

b = 0.27*2;   % wing span [m]
Ct = 0.16;   % wing tip chord [m]
Cr = 0.2;   % wing root chord [m]
c = @(y) Cr+(Ct-Cr)/(b/2)*y;    % function for chord at y
Sw = 2*integral(c,0,b/2);  % wing area [m^2]
lambda = Ct / Cr;   % tip, root chord ratio (=taper ratio)
A = (b ^ 2) / Sw;   % aspect ratio
ALE = atan((Cr-Ct)/b/2);    % leading edge sweep angle [rad]
epsilonr = 0;   % root 뒷틀림각
epsilont = 0;   % tip 뒷틀림각
diangle = 0;    % dihedral angle [rad]

%% center point of chord sweep angle: 모든 시위 중앙을 지나는 선의 뒷젖침각

centerC = Cr/2 - tan(ALE)*b/2 - (Ct / 2);    % line along center point of chord [m]
Ac2 = atan(centerC / (b / 2));  % center point of chord sweep angle [rad]

%% wing lift effectiveness: 양력선 기울기

Minf = 0.058;   % mach 0.05, 17m/s 
Beta = sqrt(1 - (Minf ^ 2));    % frandtl subsonic compressivity coefficient
cla = 0.111*180/pi;  % section lift effectiveness [/rad]     % 양력계수 기울기 (선형화 필요!!!)
kappa = cla / (2 * pi); % airfoil camber degree dimentionless coefficient [/rad]
CLa = (2 * pi * A) / (2 + sqrt(((((A ^ 2) * (Beta ^ 2)) / (kappa ^ 2)) * (1 + ((tan(Ac2) ^ 2) / Beta ^ 2)) + 4)));    % wing lift effectiveness [/rad]

%% wing zero-lift angle of attack: 날개에서 양력이 발생하지 않는 특정 받음각

a_0 =  0;    % 2D airfoil의 양력이 발생하지 않는 받음각 [rad]     %%%%%% CFD로 구해야함 (대칭형 airfoil 경우 0)
epsilon = @(y) epsilonr+(epsilont-epsilonr)/(b/2)*y;    % y에서 뒤틀림 각
fun_a_0htail = @(y) (a_0 - epsilon(y)).*c(y);    % a_0wing 계산을 위해 지정한 함수
a_0htail = 2/Sw*integral(fun_a_0htail,0,b/2);     % 3d에서 양력이 발생하지 않는 받음각 [rad]

%% 평균 공력 시위

x_LE = @(y) tan(ALE).*y;    % LE x 좌표 [m]
fun_c_bar = @(y) c(y).^2;
c_bar = 2/Sw*integral(fun_c_bar,0,b/2);     % 평균 공력시위 [m]
fun_X_LE_MAC = @(y) x_LE(y).*c(y);
X_LE_MAC = 2/Sw*integral(fun_X_LE_MAC,0,b/2);    % LE ~ MAC 까지 거리 [m]
X_AC_htail = X_LE_MAC+0.25*c_bar ;    % 공력중심 x좌표 [m]
fun_Y_MAC = @(y) y.*c(y);
Y_MAC = 2/Sw*integral(fun_Y_MAC,0,b/2); % 공력중심 y좌표 [m]

%% 날개 피칭 모멘트 계수

c_m_ac = 0;     %%%%%% 공력 중심에서의 피칭 모멘트 (양력이 0일때 모멘트 계수)  %%%% 대칭형 에어포일에서는 0
fun_C_M_AC1 = @(y) c_m_ac.*c(y).^2;
x_ac = @(y) x_LE(y)+0.25.*c(y);
fun_C_M_AC2 = @(y) cla.*(a_0htail+epsilon(y)-a_0).*(x_ac(y)-X_AC_htail).*c(y);
C_M_AC = 2/Sw/c_bar*(integral(fun_C_M_AC1,0,b/2)-integral(fun_C_M_AC2,0,b/2));  % wing pitching moment codfficient 수평꼬리날개 피칭 모멘트 계수

%% downwash 다운워시 (주익 받음각 : a_wing 꼬리날개 받음각 : a_htail)

a_wing = 2*pi/180;     % 주익 받음각 [rad]
a_0wing =  -0.10786;    % 주익 양력 0이 되는 받음각 [rad]
CLa_wing = 4.84003;       % 주익 CLa [rad]
CLa_0wing = CLa_wing*(0-a_0wing); % 받음각 0일때 양력계수
A_eff = 8.6957;          % 주익 가로세로비
l_2= 0.7-0.26;         % 주익 끝 ~ 꼬리날개 앞  [m]
b_w = 2;         % 주익 길이

e0 = 2*CLa_0wing/pi/A_eff;
de_over_da_inf = 2*CLa_wing/pi/A_eff;     % infinity downwash gradient
tail_length_in_semispans = 2*l_2/b_w;
downwash_a_w = 0.43;    % 수평 꼬리날개 다운워시 각 계수
di_h = 0;   % 수평 꼬리날개 incidence angle   [rad]
downwash = downwash_a_w*(a_wing-a_0wing);   % 수평 꼬리날개 다운워시 각 [rad]
a_htail = a_wing-downwash+di_h;     % 수평 꼬리날개 받음각   [rad]

%% 롤링 모멘트 계수 (주익 받음각 : a_wing일 때, side slip angle : B일때)

B = 2*pi/180;  % side slip 옆 미끄러짐각 [rad]

C_L_roll_ALE = CLa*(a_htail - a_0htail)*Y_MAC/b*(-sin(2*ALE)*sin(2*B));    % 뒷젖침각과 옆 미끄러짐각으로 생기는 롤링 모멘트 계수
C_L_roll_diangle = -diangle*CLa*B*Y_MAC/b;   % 상반각과 옆 미끄러워짐각으로 생기는 롤링 모멘트 계수
C_L_roll = C_L_roll_ALE+C_L_roll_diangle;    % 수평꼬리날개의 롤링 모멘트 [단위없음]

%% 항력계수 (주익 받음각 a_wing일 때)

t_c = 0.10;  % thickness ratio of the airfoil section
S_wet = 2*Sw;  % 공기에 노출된 날개 면적 (2*Sw와 비슷)      % 모델링 이후 면적 값 입력 가능
p = 1.225;  % density of air [kg/m^3]
mu = 1.54*10^(-5); % 공기 점성계수 [m^2/s]
V = 17;     % 속도
l = c(Y_MAC);   % 특성길이
R_l = p*V*l/mu; % Reynolds number

C_f = 0.005575; % skin friction coefficient 표면 마찰계수    (Figure 5.31) (aerotoolbox사이트 계산기 있음)
C_D_P = C_f*(1+2*(t_c)+100*(t_c)^4)*S_wet/Sw;   % 마찰 항력계수
C_L = CLa*(a_htail-a_0htail);     % 양력계수
e = 0.93;       % osweld's efficiency factor    (Figure osweld's effiency factor)
C_D_I = C_L^2/pi/A/e;     % 유도 항력계수
C_D = C_D_P+C_D_I;        % 항력계수 [/rad]

%% Elevator lift effectiveness : 엘리베이터 양력 효과

b_ai = b/2*0.1;  % inner coordinates of aileron [m]
b_ao = b/2*0.9;  % outer coordinates of aileron [m]
Ce_C = 0.3;  % ratio of aileron chord to wing chord

cl_delta_theory = 4.49; % Theoretical section-lift effectiveness [/rad]  (Figure 5.14)
cla_theory = 2*pi;      % theory lift effecriveness [/rad]
cla_over_cla_theory = cla/cla_theory;   % ratio of cla
cl_over_cl_theory = 1; % Empirical section-lift effectiveness  (Figure 5.15)

C_l_delta_E = 1/Beta*(cl_over_cl_theory)*cl_delta_theory; 
fun_C_L_delta_E = @(y) C_l_delta_E.*c(y);
C_L_delta_E = 2/Sw*integral(fun_C_L_delta_E,b_ai,b_ao);   % 엘리베이터 양력 효과

%% Elevator pitching moment effectiveness : 엘리베이터 피칭 모멘트 효과

C_m_del = -2*sqrt(Ce_C*(1-Ce_C)^3);     % 2d 엘리베이터 피칭 모멘트 효과 [/rad]
fun1_C_M_del_E = @(y) C_m_del.*c(y).^2;
fun2_C_M_del_E = @(y) C_l_delta_E.*(x_ac(y)-X_AC_htail).*c(y);
C_M_del_E = 2/Sw/c_bar*(integral(fun1_C_M_del_E,b_ai,b_ao)-integral(fun2_C_M_del_E,b_ai,b_ao));     % 3d 엘리베이터 피칭 모멘트 효과 [/rad]

%% Elevator drag effectiveness : 엘리베이터 항력 효과

del_E = 20*pi/180;     % 엘리베이터 변화 각 [rad]
del_c_d = 0.02;     % delta 각에 따른 항력 계수 (Figure 5.17)
fun_C_D_delta_E = @(y) del_c_d/del_E.*c(y);
C_D_delta_E = 2/Sw*integral(fun_C_D_delta_E,b_ai,b_ao);  % 엘리베이터 항력 효과

