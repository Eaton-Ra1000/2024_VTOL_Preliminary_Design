% rudder design

clear all
close all
clc


fun = @fun1;
x = fsolve(fun,[0 0]);
sigma = x(1)*180/pi % [degree]
dr = x(2)*180/pi    % [degree]

function F = fun1(x)
br = 0.27*0.9;
bv = 0.27;
Cr_C = 0.2;     %%%%% 초기 설정값
Tr = -6.624.*(Cr_C).^4+12.07.*(Cr_C).^3-8.292.*(Cr_C).^2+3.295.*(Cr_C)+0.004942

p = 1.225;
V_f = 10;       %%% 초기 설정값
V_w = 8;        %%% 초기 설정값
V_t = sqrt(V_f^2+V_w^2);
beta = atan(V_w/V_f);
S = 0.46;
Ss = 0.0486;
b = 2;
Cnbeta = 0.1139;
dc = 0.35-0.064;
Cdy = 0.7; % typically 0.55 ~ 0.8
Clav = 2.1;
Cndr = -Clav*0.0285*Tr*br/bv ;
Cybeta = -0.7*Clav*0.9*(Ss/S);
Cydr = Clav*0.9*Tr*br/bv*Ss/S;
Fw = p/2*V_w^2*Ss*Cdy;

F(1) = p/2*V_t^2*S*b*(Cnbeta*(beta-x(1))+Cndr*x(2))+Fw*dc*cos(x(1));
F(2) = p/2*V_w^2*Ss*Cdy-p/2*V_t^2*S*(Cybeta*(beta-x(1))+Cydr*x(2));
end