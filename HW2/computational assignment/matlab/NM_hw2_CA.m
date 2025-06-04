clc; clear; close all;

%% Influenc of N ==========================================================
% =========================================================================
% =========================================================================

parameters
% Ns = [5, 10, 15, 20, 40, 60, 100];
Ns = [3, 5, 7, 10];

results_vec = {};
r_vec  = {};
for Ns_index = 1:length(Ns)
    N       = Ns(Ns_index);
    h       = (r_max - r_min) / (N+1);
    R       = 0.45;
    delta_t = R * h^2;
    r       = r_min+h*(0:1:N+1);

    results_vec{Ns_index,1} = solver(r_min, r, K, h, delta_t, t_start, t_end, N, epsilon, max_iteration);
end
%%
load("effect_of_N_to_60.mat")

fig1 = figure('Name', '1','Position', [50, 250, 900, 600]);
size = 15;

colors = cool(length(results_vec(:,1)))*0.9;
for i = length(results_vec(:,1))
    plot(results_vec{i,1}.r, results_vec{i,1}.Ts(end,:), 'LineWidth', 2, 'Color', colors(i,:))
end
xlabel('r$\displaystyle\left[\frac{m}{sec}\right]$', 'FontSize',size, 'Interpreter','latex')
ylabel('T$\left[K\right]$', 'FontSize',size, 'Interpreter','latex')
title('Temperature as a Function of The Radius', 'FontSize',size, 'Interpreter','latex')
subtitle(sprintf('$N=%g$ $||$ $t=%f[s]$ $||$ $\\varepsilon=%g$ $||$ R=%g', results_vec{end,1}.N, results_vec{end,1}.t_vec(end), results_vec{end,1}.epsilon, results_vec{end,1}.R), 'FontSize', size, 'Interpreter','latex')
grid on
grid minor
box on
% exportgraphics(fig1, 'images/T over r.png','Resolution',400);

% /////////////////////////////////////////////////////////////////////////
% /////////////////////////////////////////////////////////////////////////

fig2 = figure('Name', '2','Position', [150, 250, 900, 600]);
size = 15;

colors = cool(length(results_vec(:,1)))*0.9;
lg = {};
for i = 1:length(results_vec(:,1))
    semilogy(results_vec{i,1}.r, results_vec{i,1}.Ts(end,:), 'LineWidth', 2, 'Color', colors(i,:))
    hold on;
    lg{end+1} = sprintf('N=%g, $\\Delta t=%f$', results_vec{i,1}.N, results_vec{i,1}.delta_t);
end
xlabel('r$\displaystyle\left[\frac{m}{sec}\right]$', 'FontSize',size, 'Interpreter','latex')
ylabel('T$\left[K\right]$', 'FontSize',size, 'Interpreter','latex')
title('Temperature as a Function of The Radius For Different N', 'FontSize',size, 'Interpreter','latex')
subtitle(sprintf('$t=%f[s]$ $||$ $\\varepsilon=%g$ $||$ R=%g', results_vec{end,1}.t_vec(end), results_vec{end,1}.epsilon, results_vec{end,1}.R), 'FontSize', size, 'Interpreter','latex')
legend(lg, 'FontSize',size-2, 'Location','eastoutside','Interpreter','latex')
grid on
grid minor
box on
% exportgraphics(fig1, 'images/T over r.png','Resolution',400); exportgraphics(fig2, 'images/Influenc of N.png','Resolution',400);

% /////////////////////////////////////////////////////////////////////////
% /////////////////////////////////////////////////////////////////////////

fig3 = figure('Name', '3','Position', [250, 250, 900, 600]);
rms_vec = [];
Ns      = [];
for i = 1:length(results_vec(:,1))-1
    rms_vec(i) = RMS(results_vec{i,1}.r, results_vec{i,1}.Ts(end,:), results_vec{i+1,1}.r, results_vec{i+1,1}.Ts(end,:));
    Ns(i) = results_vec{i,1}.N;
end

plot(Ns, rms_vec, 'LineWidth', 2, 'Color', colors(end,:))
    
xlabel('N', 'FontSize',size, 'Interpreter','latex')
ylabel('rms$\left[K\right]$', 'FontSize',size, 'Interpreter','latex')
title('The RMS between two solution as a function of N', 'FontSize',size, 'Interpreter','latex')
subtitle(sprintf('$t=%f[s]$ $||$ $\\varepsilon=%g$ $||$ R=%g', results_vec{end,1}.t_vec(end), results_vec{end,1}.epsilon, results_vec{end,1}.R), 'FontSize', size, 'Interpreter','latex')
grid on
grid minor
box on
% exportgraphics(fig1, 'images/T over r.png','Resolution',400); exportgraphics(fig2, 'images/Influenc of N.png','Resolution',400); exportgraphics(fig3, 'images/Influenc of N - error.png','Resolution',400);

%% Influenc of epsilon ====================================================
% =========================================================================
% =========================================================================

parameters

epsilon_vec = logspace(-1, -10, 10);
% epsilon_vec = logspace(-1, -6, 6);

results_vec = {};
r_vec  = {};
for eps_index = 1:length(epsilon_vec)
    epsilon = epsilon_vec(eps_index);
    results_vec{eps_index,1} = solver(r_min, r, K, h, delta_t, t_start, t_end, N, epsilon, max_iteration);
end
%%
load("effect_of_epsilon_-1_to_-10.mat")

fig4 = figure('Name', '4','Position', [350, 250, 900, 600]);
size = 15;

colors = cool(length(results_vec(:,1)))*0.9;
lg = {};
for i = 1:length(results_vec(:,1))
    semilogy(results_vec{i,1}.r, results_vec{i,1}.Ts(end,:), 'LineWidth', 2, 'Color', colors(i,:))
    hold on;
    lg{end+1} = sprintf('$\\varepsilon=%g$', results_vec{i,1}.epsilon);
end
xlabel('r$\displaystyle\left[\frac{m}{sec}\right]$', 'FontSize',size, 'Interpreter','latex')
ylabel('T$\left[K\right]$', 'FontSize',size, 'Interpreter','latex')
title('Temperature as a Function of The Radius For Different $\varepsilon$', 'FontSize',size, 'Interpreter','latex')
subtitle(sprintf('$t=%f[s]$ $||$ $N=%g$ $||$ R=%g', results_vec{end,1}.t_vec(end), results_vec{end,1}.N, results_vec{end,1}.R), 'FontSize', size, 'Interpreter','latex')
legend(lg, 'FontSize',size-2, 'Location','eastoutside','Interpreter','latex')
grid on
grid minor
box on
% exportgraphics(fig4, 'images/Influenc of epsilon.png','Resolution',400);

fig5 = figure('Name', '5','Position', [450, 250, 900, 600]);
rms_vec = [];
eps_vec = [];
for i = 1:length(results_vec(:,1))-1
    rms_vec(i) = RMS(results_vec{i,1}.r, results_vec{i,1}.Ts(end,:), results_vec{i+1,1}.r, results_vec{i+1,1}.Ts(end,:));
    eps_vec(i) = results_vec{i,1}.epsilon;
end

loglog(eps_vec, rms_vec, 'LineWidth', 2, 'Color', colors(end,:))
    
xlabel('$\varepsilon$', 'FontSize',size, 'Interpreter','latex')
ylabel('rms$\left[K\right]$', 'FontSize',size, 'Interpreter','latex')
title('The RMS between two solution as a function of $\varepsilon$', 'FontSize',size, 'Interpreter','latex')
subtitle(sprintf('$t=%f[s]$ $||$ $N=%g$ $||$ R=%g', results_vec{end,1}.t_vec(end), results_vec{end,1}.N, results_vec{end,1}.R), 'FontSize', size, 'Interpreter','latex')
grid on
grid minor
box on
% exportgraphics(fig4, 'images/Influenc of epsilon.png','Resolution',400); exportgraphics(fig5, 'images/Influenc of epsilon - error.png','Resolution',400);

%% Influenc of R ==========================================================
% =========================================================================
% =========================================================================

parameters

% Rs = linspace(0.05,0.49,10);  
Rs = linspace(0.45, 0.3, 4);

results_vec = {};
r_vec  = {};
for R_index = 1:length(Rs)
    R = Rs(R_index);
    delta_t = R * h^2;
    results_vec{R_index,1} = solver(r_min, r, K, h, delta_t, t_start, t_end, N, epsilon, max_iteration);
end
%%
load("effect_of_R_0.05_to_0.49.mat")

fig6 = figure('Name', '6','Position', [550, 250, 900, 600]);
size = 15;

colors = cool(length(results_vec(:,1)))*0.9;
lg = {};
for i = 1:length(results_vec(:,1))
    semilogy(results_vec{i,1}.r, results_vec{i,1}.Ts(end,:), 'LineWidth', 2, 'Color', colors(i,:))
    hold on;
    lg{end+1} = sprintf('$R=%g$', results_vec{i,1}.R);
end
xlabel('r$\displaystyle\left[\frac{m}{sec}\right]$', 'FontSize',size, 'Interpreter','latex')
ylabel('T$\left[K\right]$', 'FontSize',size, 'Interpreter','latex')
title('Temperature as a Function of The Radius For Different R', 'FontSize',size, 'Interpreter','latex')
subtitle(sprintf('$t=%f[s]$ $||$ $N=%g$ $||$ $\\varepsilon=%g$', results_vec{end,1}.t_vec(end), results_vec{end,1}.N, results_vec{end,1}.epsilon), 'FontSize', size, 'Interpreter','latex')
legend(lg, 'FontSize',size-2, 'Location','eastoutside','Interpreter','latex')
grid on
grid minor
box on
% exportgraphics(fig6, 'images/Influenc of R.png','Resolution',400);

fig7 = figure('Name', '7','Position', [650, 250, 900, 600]);
rms_vec = [];
R_vec = [];
for i = 1:length(results_vec(:,1))-1
    rms_vec(i) = RMS(results_vec{i,1}.r, results_vec{i,1}.Ts(end,:), results_vec{i+1,1}.r, results_vec{i+1,1}.Ts(end,:));
    R_vec(i) = results_vec{i,1}.R;
end

semilogy(R_vec, rms_vec, 'LineWidth', 2, 'Color', colors(end,:))
    
xlabel('$R$', 'FontSize',size, 'Interpreter','latex')
ylabel('rms$\left[K\right]$', 'FontSize',size, 'Interpreter','latex')
title('The RMS between two solution as a function of R', 'FontSize',size, 'Interpreter','latex')
subtitle(sprintf('$t=%f[s]$ $||$ $N=%g$ $||$ $\\varepsilon=%g$', results_vec{end,1}.t_vec(end), results_vec{end,1}.N, results_vec{end,1}.epsilon), 'FontSize', size, 'Interpreter','latex')
grid on
grid minor
box on
% exportgraphics(fig6, 'images/Influenc of R.png','Resolution',400); exportgraphics(fig7, 'images/Influenc of R - error.png','Resolution',400);

%% Integrate ==============================================================
% =========================================================================
% =========================================================================

parameters
result = solver(r_min, r, K, h, delta_t, t_start, t_end, N, epsilon, max_iteration);

% integrate N=20
sum_N20 = 0;
for i = [0:result.N+1-1]+1
    sum_N20 = sum_N20 + (result.Ts(end,i)+result.Ts(end,i+1))*(result.r(i)+result.r(i+1))/2;
end
sum_N20 = sum_N20 * alpha * 0.5 * result.h

% integrate N=40
load("effect_of_N_to_60.mat")
result = results_vec{end-1};
sum_N40 = 0;
for i = [0:result.N+1-1]+1

    sum_N40 = sum_N40 + (result.Ts(end,i)+result.Ts(end,i+1))*(result.r(i)+result.r(i+1))/2;
end
sum_N40 = sum_N40 * alpha * 0.5 * result.h

% integrate N=60
load("effect_of_N_to_60.mat")
result = results_vec{end};
sum_N60 = 0;
for i = [0:result.N+1-1]+1

    sum_N60 = sum_N60 + (result.Ts(end,i)+result.Ts(end,i+1))*(result.r(i)+result.r(i+1))/2;
end
sum_N60 = sum_N60 * alpha * 0.5 * result.h




%% Functions ==============================================================
% =========================================================================
% =========================================================================

function T = init_fuild(r_min, h, N)
    i = 0:1:N+1;
    r = r_min+h * i;
    T = 200 * (r - 0.5);
end

function T = set_BC(T, t)
    T(1)   = t;
    T(end) = 100 + 40 * t;
end

function T_next = step_space_jacobi(T_current, r, K, h, delta_t, N)
    T_next(1) = T_current(1);
    T_next(N+1+1) = T_current(N+1+1);    
    for i = [1:N]+1
        T_next(i) = T_current(i) + 4 * K * delta_t * ((T_current(i+1) - 2 * T_current(i) + T_current(i-1)) / (h^2) + 1 / r(i) * (T_current(i+1) - T_current(i-1)) / (2*h));
        % T_next(i) = T_current(i) + 4 * K * delta_t * ((T_current(i+1) - 2 * T_current(i) + T_next(i-1)) / (h^2) + 1 / r(i) * (T_current(i+1) - T_next(i-1)) / (2*h));
    end
end

function converged = check_convergence(T_next, T_current, epsilon, N)
    converged = true;
    for i = [1:N]+1
        if abs(T_next(i) - T_current(i)) > epsilon
            converged = false;
            return
        end
    end
end

function T_next_t = solve_for_specific_t_jacobi(T_current_t, r, K, h, delta_t, N, epsilon, max_iteration)
    T_current = T_current_t;
    for n = 1:max_iteration
        T_next = step_space_jacobi(T_current, r, K, h, delta_t, N);
        if check_convergence(T_next, T_current, epsilon, N)
            break
        end
        T_current = T_next;
    end
    T_next_t = T_next;
end

function results = solver(r_min, r, K, h, delta_t, t_start, t_end, N, epsilon, max_iteration)
    T_current_t = init_fuild(r_min, h, N);
    t = t_start;
    Ts(1, :) = T_current_t;
    
    while t <= t_end
        fprintf('N: %d || R: %g || epsilon: %g || delta_t: %g || t: %6.4f/%g\n', N, delta_t / h^2, epsilon, delta_t, t, t_end)
        T_current_t = set_BC(T_current_t, t);
        T_next_t = solve_for_specific_t_jacobi(T_current_t, r, K, h, delta_t, N, epsilon, max_iteration);
        Ts(end+1, :) = T_next_t;
        
        t = t + delta_t;
        T_current_t = T_next_t;
    end
    results.Ts      = Ts;
    results.N       = N;
    results.r       = r;
    results.h       = h;
    results.delta_t = delta_t;
    results.epsilon = epsilon;
    results.R       = delta_t / h^2;
    results.t_vec   = delta_t*(0:1:t_end/delta_t);
end

function y = interpulate(x, x_vec, y_vec)
% This function assume an ordered x_vec from low to high
    if x > x_vec(end) || x < x_vec(1)
        fprintf('value out of bounds\n');
        return;
    end

    index_i = 0;
    for i = 1:length(x_vec)
        if x <= x_vec(i)
            index_i = i;
            break;
        end
    end
    
    if x == x_vec(index_i)
        y = y_vec(index_i);
        return;
    end

    m = (y_vec(index_i+1) - y_vec(index_i)) / (x_vec(index_i+1) - x_vec(index_i));
    b = y_vec(index_i) - x_vec(index_i) * m;

    y = m * x + b;
end

function rms = RMS(x1, y1, x2, y2)
% This functions calcutlates the RMS value for 1 compared to 2
    y2_len_y1 = [];
    for i = 1:length(x1)
        y2_len_y1(i) = interpulate(x1(i), x2, y2);
    end

    rms = sqrt(sum((y1-y2_len_y1).^2)/length(y1));
end


