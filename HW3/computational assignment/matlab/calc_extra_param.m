AB = 12; BC = 12; CD = 2; DE = 6 ;
EF = 6 ; FG = 6 ; GH = 4; HA = 12;

ni = lcm(BC, BC-DE) * factor;
nj = lcm(lcm(CD, AB), lcm(CD+EF, AB)) * factor;
delta_x = (x_max - x_min) / (ni);
delta_y = (y_max - y_min) / (nj);
x_vec = x_min + (0:ni)*delta_x;
y_vec = y_min + (0:nj)*delta_y;

i_AB = 0;
j_AB = 0:(AB / delta_y);

i_BC = 0:(BC / delta_x);
j_BC = 0;

i_CD = BC / delta_x;
j_CD = 0:(CD / delta_y);

i_DE = ((BC - DE) / delta_x):(BC / delta_x);
j_DE = CD / delta_y;

i_EF = (BC - DE) / delta_x;
j_EF = (CD / delta_y):((CD + EF) / delta_y);

i_FG = ((BC - DE) / delta_x):(BC / delta_x);
j_FG = (CD + EF) / delta_y;

i_GH = BC / delta_x;
j_GH = ((CD + EF) / delta_y):((CD + EF + GH) / delta_y);

i_HA = i_BC;
j_HA = (AB / delta_y);

for j = j_AB+1
    x_mat(j, :) = x_vec;
end
for i = i_BC+1
    y_mat(:, i) = y_vec;
end