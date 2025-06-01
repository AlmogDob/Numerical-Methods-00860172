N             = 20;
epsilon       = 1e-5;
max_iteration = 1e6;

r_max   = 1;
r_min   = 0.5;
t_start = 0;
t_end   = 10;
K       = 0.1;
alpha   = 10.7;

h       = (r_max - r_min) / (N+1);
R       = 0.49;
% R       = delta_t / h^2;
delta_t = R * h^2;
r       = r_min+h*(0:1:N+1);

