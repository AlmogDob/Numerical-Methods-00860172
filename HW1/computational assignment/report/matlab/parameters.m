
N       = 101;
epsilon = 1e-12;

lambda     = 0.1;
T_v        = 0.203;
T_u        = 0.152;
alpha_beta = 0.0234;
infinity   = 30;
zeta_max   = infinity;
zeta_min   = -infinity;
h          = (zeta_max - zeta_min) / (N+1);
T_start    = zeta_min * (T_v - T_u);
T_end      = zeta_max * (T_v - T_u + alpha_beta);
md_start   = 1;
md_end     = 0;
