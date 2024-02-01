clear all
close all
clc

%% vertical tail wing airfoil is NACA0010

%% 3-D wing condition and definition

b = 0.27;   % vertical tail wing span [m]
Ct = 0.16;   % vertical tail wing tip chord [m]
Cr = 0.2;   % vertical tail wing root chord [m]
c = @(z) Cr+(Ct-Cr)/(b)*z;    % function for chord at z
Sv = integral(c,0,b);  % rudder area [m^2]
lambda = Ct / Cr;   % tip, root chord ratio (=taper ratio)
A = (b ^ 2) / Sv;   % aspect ratio
ALE = atan((Cr-Ct)/b);    % leading edge sweep angle [rad]
epsilonr = 0;   % root 뒷틀림각
epsilont = 0;   % tip 뒷틀림각
diangle = 0;    % dihedral angle [rad]

%% center point of chord sweep angle: 모든 시위 중앙을 지나는 선의 뒷젖침각

centerC = Cr/2 - tan(ALE)*b - (Ct / 2);    % line along center point of chord [m]
Ac2 = atan(centerC / (b));  % center point of chord sweep angle [rad]

%% sideforce effectiveness : 횡력선 기울기

Minf = 0.058;   % mach 0.05, 17m/s 
Beta = sqrt(1 - (Minf ^ 2));    % frandtl subsonic compressivity coefficient
csa = 0.111*180/pi;  % section sideforce effectiveness [/rad]     % 횡력계수 기울기 (선형화 필요!!!)
kappa = csa / (2 * pi); % airfoil camber degree dimentionless coefficient [/rad]
C_S_B = (2 * pi * A) / (2 + sqrt(((((A ^ 2) * (Beta ^ 2)) / (kappa ^ 2)) * (1 + ((tan(Ac2) ^ 2) / Beta ^ 2)) + 4)));    % vertical tail wing sideforce effectiveness [/rad]

%% vertical tail wing zero-sideforce angle of attack: 날개에서 횡력이 발생하지 않는 특정 받음각

B_0 =  0;    % 2D airfoil의 횡력이 발생하지 않는 받음각 [rad]     %%%%%% CFD로 구해야함 (대칭형 airfoil 경우 0)
epsilon = @(z) epsilonr+(epsilont-epsilonr)/(b)*z;    % z에서 뒤틀림 각
fun_B_0vtail = @(z) (B_0 - epsilon(z)).*c(z);    % B_0tail 계산을 위해 지정한 함수
B_0vtail = 1/Sv*integral(fun_B_0vtail,0,b);     % 3d에서 횡력이 발생하지 않는 받음각 [rad]

%% 평균 공력 시위

x_LE = @(z) tan(ALE).*z;    % LE x 좌표 [m]
fun_c_bar = @(z) c(z).^2;
c_bar = 1/Sv*integral(fun_c_bar,0,b);     % 평균 공력시위 [m]
fun_X_LE_MAC = @(z) x_LE(z).*c(z);
X_LE_MAC = 1/Sv*integral(fun_X_LE_MAC,0,b);    % LE ~ MAC 까지 거리 [m]
X_AC_vtail = X_LE_MAC+0.25*c_bar ;    % 공력중심 x좌표 [m]
fun_Z_MAC = @(z) z.*c(z);
Z_MAC = 1/Sv*integral(fun_Z_MAC,0,b); % 공력중심 z좌표 [m]
x_ac = @(z) x_LE(z)+0.25.*c(z);

%% 항력계수 (받음각 B_vtail일 때)

B_vtail = 2*pi/180;     % 받음각 [rad]
t_c = 0.10;  % thickness ratio of the airfoil section
S_wet = 2*Sv;  % 공기에 노출된 날개 면적 (2*Sw와 비슷)      % 모델링 이후 면적 값 입력 가능
p = 1.225;  % density of air [kg/m^3]
mu = 1.54*10^(-5); % 공기 점성계수 [m^2/s]
V_inf = 17;     % 속도 [m/s]
l = c(Z_MAC);   % 특성길이
R_l = p*V_inf*l/mu; % Reynolds number

C_f = 0.005575; % skin friction coefficient 표면 마찰계수    (Figure 5.31) (aerotoolbox사이트 계산기 있음)
C_D_P = C_f*(1+2*(t_c)+100*(t_c)^4)*S_wet/Sv;   % 마찰 항력계수
C_S = C_S_B*(B_vtail-B_0vtail);     % 횡력계수
e = 0.93;       % osweld's efficiency factor    (Figure osweld's effiency factor)
C_D_I = C_S^2/pi/A/e;     % 유도 항력계수
C_D = C_D_P+C_D_I;        % 항력계수 [/rad]

%% Rudder sideforce effectiveness : 러더 횡력 효과

b_ai = b*0.1;  % inner coordinates of rudder [m]
b_ao = b*0.9;  % outer coordinates of rudder [m]
Ce_C = 0.3;  % ratio of rudder chord to vertical tail wing chord

cs_delta_theory = 4.49; % Theoretical section-sideforce effectiveness [/rad]  (Figure 5.14)
csa_theory = 2*pi;      % theory sideforce effecriveness [/rad]
csa_over_csa_theory = csa/csa_theory;   % ratio of cla
cs_over_cs_theory = 1; % Empirical section-sideforce effectiveness  (Figure 5.15)

C_s_delta_R = 1/Beta*(cs_over_cs_theory)*cs_delta_theory;
fun_C_S_delta_R = @(z) C_s_delta_R.*c(z);
C_S_delta_R = 1/Sv*integral(fun_C_S_delta_R,b_ai,b_ao);   % 러더 횡력 효과 [/rad]

%% Rudder yawing moment effectiveness : 러더 요잉 모멘트 효과

C_m_del = -1*sqrt(Ce_C*(1-Ce_C)^3);     % 2d 러더 요잉 모멘트 효과
fun1_C_M_del_R = @(z) C_m_del.*c(z).^2;
fun2_C_M_del_R = @(z) C_s_delta_R.*(x_ac(z)-X_AC_vtail).*c(z);
C_M_del_R = 1/Sv/c_bar*(integral(fun1_C_M_del_R,b_ai,b_ao)-integral(fun2_C_M_del_R,b_ai,b_ao));     % 3d 러더 요잉 모멘트 효과 [/rad]

%% Rudder drag effectiveness : 러더 항력 효과

del_R = 20*pi/180;     % 러더 변화 각 [rad]
del_c_d = 0.02;     % delta 각에 따른 항력 계수 (Figure 5.17)
fun_C_D_delta_R = @(z) del_c_d/del_R.*c(z).*z;
C_D_delta_R = 1/Sv*integral(fun_C_D_delta_R,b_ai,b_ao);  % 엘리베이터 항력 효과 [/rad]

