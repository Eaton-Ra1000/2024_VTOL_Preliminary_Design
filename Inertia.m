%% Mass moments of inerita (설계단계에서 관성모멘트 가정)
b = 2;         % [m]
W = 6*9.81;    % [kg*m/s^2]
g = 9.81;      % [m/s^2]
L = 0.7;     % [m]
R_x = 0.25;
R_y = 0.38;
R_z = 0.39;

I_xx = b^2*W*R_x^2/4/g;  % [kg*m^2]
I_yy = L^2*W*R_y^2/4/g;
I_zz = ((b+L)/2)^2*W*R_z^2/4/g;