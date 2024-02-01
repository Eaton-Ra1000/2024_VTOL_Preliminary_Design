clear all
close all
clc
format shortE

%% basic factor

L = 1000;                   % quad size [mm]
La = (L / sqrt(2)) / 2;     % quad 1/2 horizontal length [mm]
Lb = L / sqrt(2);           % quad longitudinal length [mm]
x = -30;                      % quad c.g margin [mm]
percentMAC = 0.25;          % %MAC

d = 13;                     % quad prop size [inch]
r = (d * 25.4) / 2;         % radius of quad prop [mm]
P = 5.5;                    % quad prop pitch [inch]
dp = 12;                    % puller prop size [inch]
rp = (dp * 25.4) / 2;       % radius of puller prop [mm]
Pp = 7;                     % puller prop pitch [inch]

Ws = 1800;                  % wing span [mm]
Ls = Ws / 2;                % half span [mm]
Wt = 540;                   % tail span [mm]
Ts = Wt / 2;                % tail half span [mm]

Cr = 260;                   % wing root chord [mm]
Ct = 220;                   % wing tip chord [mm]
Tr = 200;                   % tail root chord [mm]
Tt = 160;                   % tail tip chord [mm]
W = 120;                    % fuselage width [mm]
D = 700;                    % distance wing leading edge to tail [mm]

m = 6;                      % aircraft mass [kg]
g = 9.81;                   % standard gravity [m/s^2]
rho = 1.225;                % air density [kg/m^3]
Cl = 1.2;                   % lift coefficient
Cd = 0.015;                 % drag coefficient
alpha = 5;                  % angle of attack [degree]

S = (Ct + Cr) * Ls * (10^(-6));                          % wing area(tapered) [m^2]
t = Ct / Cr;                                             % taper ratio
MAC = Cr * (2 / 3) * ((1 + t + t^2) / (1 + t));          % mean aerodynamic chord [mm]
Lac = MAC * percentMAC;                                  % aerodynamic center from leading edge [mm]
AR = Ws / MAC;                                           % aspect ratio
Cx = Cr - ((Cr - Ct) / Ls) * (Ls - La)                   % chord when wing is perfect taper [mm]
Vs = sqrt((2 * m * g) / (rho * S * Cl))                  % stall speed [m/s]
WCL = m / (S * sqrt(S))                                  % cubic wing loading [kg/m^3]
LDR = Cl / Cd                                            % lift drag ratio

tt = Tt / Tr;                                            % tail taper ratio
MACt = Tr * (2 / 3) * ((1 + tt + tt^2) / (1 + tt));      % tail mean aerodynamic chord [mm]
St = (Tt + Tr) * Ts * (10^-6);                           % tail area[mm^2]
Lact = MACt * percentMAC;                                % tail aerodynamic center [mm]
Ta = (D - Lac) + Lact;                                   % tail arm [mm]
Vbar = (St / S) * (Ta / MAC);                            % tail volume
Np = percentMAC + (percentMAC * sqrt(sqrt(AR))) * Vbar;  % neutral point [%MAC]
Lnp = Np * MAC                                           % neutral point [mm]
MAC_x = (Ws/2)*(1-(MAC-Ct)/(Cr-Ct));                     % MAC(x)

result = [];

for Lcg = Cr / 2 : -1 : Lac           % c.g from leading edge [mm]
    Ln = La - Lcg - x;                % nose overhnag length [mm]
 
    result = [result; Lcg, MAC, Lac, S];
end
 
disp('Lcg | MAC | Lac | S');
disp(result);

%% chosen cg

xcg = 40;                % chosen c.g from leading edge [mm]
cg = xcg - x;            % chosen c.g with c.g margin [mm]

%% eVTOL drawing

TailStart = Cr/2-D;
WingVertices = [Ls, Ct/2;...        % point 1
                Ls, -Ct/2;...       % point 2
                -Ls, -Ct/2;...      % point 3
                -Ls, Ct/2;...       % point 4
                0, Cr/2;...         % point 5
                0, -Cr/2;...        % point 6

                0, TailStart;...            % point 7    
                0, TailStart-Tr;...         % point 8
                Wt/2, TailStart-Tr;...      % point 9
                Wt/2, TailStart-Tr+Tt       % point 10 
                -Wt/2, TailStart-Tr;...     % point 11
                -Wt/2, TailStart-Tr+Tt;...  % point 12

                W/2, ((Cr / 2) + Ln);...    % point 13
                -W/2, ((Cr / 2) + Ln);...   % point 14
                W/2, TailStart-Tr;...       % point 15
                -W/2, TailStart-Tr;...      % point 16
                ];

WingFaces = [5, 1, 2, 6, 3, 4;
             7, 8, 9, 10, 10, 10;
             7, 8, 11, 12, 12, 12;
             13, 14, 16, 15, 15, 15];

xaxis = [1, 0];
yaxis = [0, 1];
varphi = linspace(0, 2*pi, 90);
xdots = r * cos(varphi);
ydots = r * sin(varphi);
PropDots = [xdots; ydots];
affine = [1, 0;...
          0, 0.3];

Prop1 = PropDots + [La; (La + ((Cr / 2) - cg))];
Prop2 = PropDots + [-La; (-La + ((Cr / 2) - cg))];
Prop3 = PropDots + [-La; (La + ((Cr / 2) - cg))];
Prop4 = PropDots + [La; (-La + ((Cr / 2) - cg))];
Puller = affine*PropDots + [0; (((Cr / 2) + Ln))];


TWOSTONE_1 = figure(1);
TWOSTONE_1.Position = 150 * [6.2420 1.65 5 5];
patch('Vertices', WingVertices, 'Faces', WingFaces,'EdgeColor','black',...
                                             'FaceColor','none');

hold on
prop1 = scatter(Prop1(1,:), Prop1(2,:), 'r.');
prop2 = scatter(Prop2(1,:), Prop2(2,:), 'g.');
prop3 = scatter(Prop3(1,:), Prop3(2,:), 'b.');
prop4 = scatter(Prop4(1,:), Prop4(2,:), 'c.');
puller = scatter(Puller(1,:), Puller(2,:), '.');

right_MAC = plot([MAC_x, MAC_x], [-MAC/2, MAC/2], 'Color',[0 0 0]);
left_MAC = plot([-MAC_x, -MAC_x], [-MAC/2, MAC/2], 'Color',[0 0 0]);
cent_MAC = plot([-MAC_x, MAC_x], [MAC * percentMAC, MAC * percentMAC], 'Color',[0 0 0]);

scatter(0, (Cr / 2) - cg, 'ko');
scatter(0, (Cr / 2) - Lnp, 'go');
scatter(0, MAC * percentMAC, 'mo');
scatter(0, (Cr - D - Tt + (1-MAC * percentMAC)), 'co');

legend([prop1, prop2, prop3, prop4], {'prop1', 'prop2', 'prop3', 'prop4'})
xlim([-1000, 1000])
ylim([-1000, 1000])