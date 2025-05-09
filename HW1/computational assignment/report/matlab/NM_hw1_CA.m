clc; clear; close all;

%%
% Influence of Parameters
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
parameters

fig1 = figure ('Position',[0 50 900 500]);
hold all

Ns = 1:20:201;
result = {};
lg = {};
for i = 1:length(Ns)
    N = Ns(i);
    colors = cool(length(Ns));
    h = (zeta_max - zeta_min) / (N+1);
    result{end+1} = finite_difference_method(T_start, T_end, lambda, alpha_beta, T_v, T_u, epsilon, h, N);
    plot(linspace(zeta_min, zeta_max, N+2), result{end}{end,:}, '-', 'LineWidth', 1.5, 'Color', colors(i,:))
    lg{end+1} = sprintf('N = %d', N);
end

size = 20;
title('FD - Temperature as a Function of $\zeta$ for different N','FontSize',size, 'Interpreter','latex');
subtitle(sprintf('$\\varepsilon=%g$', epsilon), 'Interpreter','latex')
ylabel('$T$ $[-]$','FontSize',size, 'Interpreter','latex')
xlabel('$\zeta$ $[-]$','FontSize',size, 'Interpreter','latex')
grid on
grid minor
box on
legend(lg,'FontSize',size-3 ,'Location','eastoutside', 'Interpreter','latex')

zoom = axes('position',[0.175 0.6 0.275 0.275]);
box on % put box around new pair of axes
hold all
for i = 1:length(Ns)
    N = Ns(i);
    colors = cool(length(Ns));
    h = (zeta_max - zeta_min) / (N+1);
    plot(linspace(zeta_min, zeta_max, N+2), result{i}{end,:}, '-', 'LineWidth', 1.5, 'Color', colors(i,:))
end
zoom.XLim = [-17.4655, -17.4605];
zoom.YLim = [-0.9958, -0.99545];
grid on
grid minor
% exportgraphics(fig1, 'images/FD - T vs zeta for diff N.png','Resolution',400);

%% 

parameters

fig2 = figure ('Position',[150 50 900 500]);
hold all

Ns = 1:20:201;
ns = [];
for i = 1:length(Ns)
    N = Ns(i);
    colors = cool(length(Ns));
    h = (zeta_max - zeta_min) / (N+1);
    T_history_FD = finite_difference_method(T_start, T_end, lambda, alpha_beta, T_v, T_u, epsilon, h, N);
    ns(i) = length(T_history_FD(:,1));
end
size = 20;
plot(Ns, ns, '-', 'LineWidth', 1.5, 'Color', cool(1)*0.8)
title('FD - Number of Iterations as a Function of Number of Elements','FontSize',size, 'Interpreter','latex');
subtitle(sprintf('$\\varepsilon=%g$', epsilon), 'Interpreter','latex')
ylabel('$n$ $[-]$','FontSize',size, 'Interpreter','latex')
xlabel('$N$ $[-]$','FontSize',size, 'Interpreter','latex')
grid on
grid minor
box on

% exportgraphics(fig2, 'images/FD - n vs N.png','Resolution',400);

%% 

parameters

fig3 = figure ('Position',[300 50 900 500]);

epss = logspace(-3,-12,10);
ns = [];
for i = 1:length(epss)
    epsilon = epss(i);
    colors = cool(length(epss));
    T_history_FD = finite_difference_method(T_start, T_end, lambda, alpha_beta, T_v, T_u, epsilon, h, N);
    ns(i) = length(T_history_FD(:,1));
end
size = 20;
semilogx(epss, ns, '-', 'LineWidth', 1.5, 'Color', cool(1)*0.8)
title('FD - Number of Iterations as a Function of epsilon','FontSize',size, 'Interpreter','latex');
subtitle(sprintf('$N=%d$', N), 'Interpreter','latex')
ylabel('$n$ $[-]$','FontSize',size, 'Interpreter','latex')
xlabel('$\varepsilon$ $[-]$','FontSize',size, 'Interpreter','latex')
grid on
grid minor
box on

% exportgraphics(fig3, 'images/FD - n vs epsilon.png','Resolution',400);

%%
parameters

fig4 = figure ('Position',[450 50 900 500]);
hold all

epss = logspace(-1,-13,13);
result = {};
lg = {};
for i = 1:length(epss)
    epsilon = epss(i);
    colors = cool(length(epss));
    result{end+1} = finite_difference_method(T_start, T_end, lambda, alpha_beta, T_v, T_u, epsilon, h, N);
    plot(linspace(zeta_min, zeta_max, N+2), result{end}{end,:}, '-', 'LineWidth', 1.5, 'Color', colors(i,:))
    lg{end+1} = sprintf('$\\varepsilon$ = %d', epsilon);
end

size = 20;
title('FD - Temperature as a Function of $\zeta$ for different $\varepsilon$','FontSize',size, 'Interpreter','latex');
subtitle(sprintf('$N=%d$', N), 'Interpreter','latex')
ylabel('$T$ $[-]$','FontSize',size, 'Interpreter','latex')
xlabel('$\zeta$ $[-]$','FontSize',size, 'Interpreter','latex')
grid on
grid minor
box on
legend(lg,'FontSize',size-3 ,'Location','eastoutside', 'Interpreter','latex')

zoom = axes('position',[0.45 0.2 0.2 0.2]);
box on % put box around new pair of axes
for i = 1:length(epss)
    epsilon = epss(i);
    colors = cool(length(epss));
    plot(linspace(zeta_min, zeta_max, N+2), result{i}{end,:}, '-', 'LineWidth', 1.5, 'Color', colors(i,:))
    hold all
end
zoom.XLim = [-0.377688, -0.37768775];
zoom.YLim = [-0.02465371, -0.02465369];
grid on
grid minor

zoom = axes('position',[0.175 0.6 0.2 0.2]);
box on % put box around new pair of axes
for i = 1:length(epss)
    epsilon = epss(i);
    colors = cool(length(epss));
    plot(linspace(zeta_min, zeta_max, N+2), result{i}{end,:}, '-', 'LineWidth', 1.5, 'Color', colors(i,:))
    hold all
end
zoom.XLim = [-1.002, -0.986];
zoom.YLim = [-0.0648, -0.0636];
grid on
grid minor
% exportgraphics(fig4, 'images/FD - T vs zeta for diff epsilon.png','Resolution',400);

%%
parameters

fig5 = figure ('Position',[0 200 900 500]);
hold all

Ns = 1:20:201;
s0 = 0;
result = {};
lg = {};
for i = 1:length(Ns)
    N = Ns(i);
    colors = cool(length(Ns));
    h = (zeta_max - zeta_min) / (N+1);
    result{end+1} = shooting_method_with_zero(T_start, T_end, s0, lambda, alpha_beta, T_v, T_u, epsilon, h, N);
    plot(linspace(zeta_min, zeta_max, N+2), result{end}{end,:}, '-', 'LineWidth', 1.5, 'Color', colors(i,:))
    lg{end+1} = sprintf('N = %d', N);
end

size = 20;
title('Shooting - Temperature as a Function of $\zeta$ for different N','FontSize',size, 'Interpreter','latex');
subtitle(sprintf('$\\varepsilon=%g$, $|$ $s_0=%f$', epsilon, s0), 'Interpreter','latex')
ylabel('$T$ $[-]$','FontSize',size, 'Interpreter','latex')
xlabel('$\zeta$ $[-]$','FontSize',size, 'Interpreter','latex')
grid on
grid minor
box on
legend(lg,'FontSize',size-3 ,'Location','eastoutside', 'Interpreter','latex')

zoom = axes('position',[0.43 0.2 0.22 0.22]);
box on % put box around new pair of axes
hold all
for i = 1:length(Ns)
    N = Ns(i);
    colors = cool(length(Ns));
    h = (zeta_max - zeta_min) / (N+1);
    plot(linspace(zeta_min, zeta_max, N+2), result{i}{end,:}, '-', 'LineWidth', 1.5, 'Color', colors(i,:))
end
zoom.XLim = [-5.39, -5.35];
zoom.YLim = [-0.3415, -0.3385];
grid on
grid minor
% exportgraphics(fig5, 'images/shooting - T vs zeta for diff N.png','Resolution',400);

%% 

parameters

fig6 = figure ('Position',[150 200 900 500]);
hold all

Ns = 1:20:201;
s0 = 0;
ns = [];
for i = 1:length(Ns)
    N = Ns(i);
    colors = cool(length(Ns));
    h = (zeta_max - zeta_min) / (N+1);
    T_history_FD = shooting_method_with_zero(T_start, T_end, s0, lambda, alpha_beta, T_v, T_u, epsilon, h, N);
    ns(i) = length(T_history_FD(:,1));
end
size = 20;
plot(Ns, ns, '-', 'LineWidth', 1.5, 'Color', cool(1)*0.8)
title('Shooting - Number of Iterations as a Function of Number of Elements','FontSize',size, 'Interpreter','latex');
subtitle(sprintf('$\\varepsilon=%g$, $|$ $s_0=%f$', epsilon, s0), 'Interpreter','latex')
ylabel('$n$ $[-]$','FontSize',size, 'Interpreter','latex')
xlabel('$N$ $[-]$','FontSize',size, 'Interpreter','latex')
grid on
grid minor
box on
ylim([0,13])

% exportgraphics(fig6, 'images/shooting - n vs N.png','Resolution',400);

%% 

parameters

fig7 = figure ('Position',[300 200 900 500]);

epss = logspace(-1,-13,13);
s0 = 0;
ns = [];
for i = 1:length(epss)
    epsilon = epss(i);
    colors = cool(length(epss));
    T_history_shooting = shooting_method_with_zero(T_start, T_end, s0, lambda, alpha_beta, T_v, T_u, epsilon, h, N);
    ns(i) = length(T_history_shooting(:,1));
end
size = 20;
semilogx(epss, ns, '-', 'LineWidth', 1.5, 'Color', cool(1)*0.8)
title('Shooting - Number of Iterations as a Function of Epsilon','FontSize',size, 'Interpreter','latex');
subtitle(sprintf('$N=%d$, $|$ $s_0=%f$', N, s0), 'Interpreter','latex')
ylabel('$n$ $[-]$','FontSize',size, 'Interpreter','latex')
xlabel('$\varepsilon$ $[-]$','FontSize',size, 'Interpreter','latex')
grid on
grid minor
box on

% exportgraphics(fig7, 'images/shooting - n vs epsilon.png','Resolution',400);

%%
parameters

fig8 = figure ('Position',[450 200 900 500]);
hold all

epss = logspace(-1,-13,13);
s0 = 0;
result = {};
lg = {};
for i = 1:length(epss)
    epsilon = epss(i);
    colors = cool(length(epss));
    result{end+1} = shooting_method_with_zero(T_start, T_end, s0, lambda, alpha_beta, T_v, T_u, epsilon, h, N);
    plot(linspace(zeta_min, zeta_max, N+2), result{end}{end,:}, '-', 'LineWidth', 1.5, 'Color', colors(i,:))
    lg{end+1} = sprintf('$\\varepsilon$ = %d', epsilon);
end

size = 20;
title('Shooting - Temperature as a Function of $\zeta$ for different $\varepsilon$','FontSize',size, 'Interpreter','latex');
subtitle(sprintf('$N=%d$, $|$ $s_0=%f$', N, s0), 'Interpreter','latex')
ylabel('$T$ $[-]$','FontSize',size, 'Interpreter','latex')
xlabel('$\zeta$ $[-]$','FontSize',size, 'Interpreter','latex')
grid on
grid minor
box on
legend(lg,'FontSize',size-3 ,'Location','eastoutside', 'Interpreter','latex')

zoom = axes('position',[0.43 0.2 0.22 0.22]);
box on % put box around new pair of axes
hold all
for i = 1:length(epss)
    epsilon = epss(i);
    colors = cool(length(epss));
    plot(linspace(zeta_min, zeta_max, N+2), result{i}{end,:}, '-', 'LineWidth', 1.5, 'Color', colors(i,:))
end
zoom.XLim = [-8e-6, -6e-6];
zoom.YLim = [-8e-7, -3e-7];
grid on
grid minor
% exportgraphics(fig8, 'images/shooting - T vs zeta for diff epsilon.png','Resolution',400);

%%
parameters

fig9 = figure ('Position',[600 200 900 500]);
hold all

s0s = -0.5:0.1:0.5;
result = {};
lg = {};
for i = 1:length(s0s)
    s0 = s0s(i);
    colors = cool(length(s0s));
    result{end+1} = shooting_method_with_zero(T_start, T_end, s0, lambda, alpha_beta, T_v, T_u, epsilon, h, N);
    plot(linspace(zeta_min, zeta_max, N+2), result{end}{end,:}, '-', 'LineWidth', 1.5, 'Color', colors(i,:))
    lg{end+1} = sprintf('$s_0$ = %7.4f', s0);
end

size = 20;
title('Shooting - Temperature as a Function of $\zeta$ for different $s_0$','FontSize',size, 'Interpreter','latex');
subtitle(sprintf('$N=%d$, $|$ $\\varepsilon=%g$', N, epsilon), 'Interpreter','latex')
ylabel('$T$ $[-]$','FontSize',size, 'Interpreter','latex')
xlabel('$\zeta$ $[-]$','FontSize',size, 'Interpreter','latex')
grid on
grid minor
box on
legend(lg,'FontSize',size-3 ,'Location','eastoutside', 'Interpreter','latex')

zoom = axes('position',[0.43 0.2 0.22 0.22]);
box on % put box around new pair of axes
hold all
for i = 1:length(s0s)
    s0 = s0s(i);
    colors = cool(length(s0s));
    plot(linspace(zeta_min, zeta_max, N+2), result{i}{end,:}, '-', 'LineWidth', 1.5, 'Color', colors(i,:))
end
zoom.XLim = [-6.905842e-6, -6.905838e-6];
zoom.YLim = [-4.5194355e-7, -4.519434e-7];
grid on
grid minor
% exportgraphics(fig9, 'images/shooting - T vs zeta for diff s0.png','Resolution',400);
% exportgraphics(fig8, 'images/shooting - T vs zeta for diff epsilon.png','Resolution',400); exportgraphics(fig7, 'images/shooting - n vs epsilon.png','Resolution',400); exportgraphics(fig6, 'images/shooting - n vs N.png','Resolution',400); exportgraphics(fig5, 'images/shooting - T vs zeta for diff N.png','Resolution',400); exportgraphics(fig4, 'images/FD - T vs zeta for diff epsilon.png','Resolution',400); exportgraphics(fig3, 'images/FD - n vs epsilon.png','Resolution',400); exportgraphics(fig2, 'images/FD - n vs N.png','Resolution',400); exportgraphics(fig1, 'images/FD - T vs zeta for diff N.png','Resolution',400);

%% 
parameters

fig10 = figure ('Position',[750 200 900 500]);

ns = [];
s0s = 0:0.0001:0.3;
lg = {};
for i = 1:length(s0s)
    s0 = s0s(i);
    colors = cool(length(s0s));
    T_history_shooting = shooting_method_with_zero(T_start, T_end, s0, lambda, alpha_beta, T_v, T_u, epsilon, h, N);
    ns(i) = length(T_history_shooting(:,1));
end
size = 20;
plot(s0s, ns, '-', 'LineWidth', 1.5, 'Color', cool(1)*0.8)
title('Shooting - Number of Iterations as a Function of $s_0$','FontSize',size, 'Interpreter','latex');
subtitle(sprintf('$N=%d$, $|$ $\\varepsilon=%g$', N, epsilon), 'Interpreter','latex')
ylabel('$n$ $[-]$','FontSize',size, 'Interpreter','latex')
xlabel('$s_0$ $[-]$','FontSize',size, 'Interpreter','latex')
grid on
grid minor
box on

% exportgraphics(fig10, 'images/shooting - n vs s0 - no negative.png','Resolution',400);

%% 
parameters

fig11 = figure ('Position',[750 300 900 500]);

ns = [];
s0s = -0.2:0.01:0.5;
lg = {};
for i = 1:length(s0s)
    s0 = s0s(i);
    colors = cool(length(s0s));
    T_history_shooting = shooting_method_with_zero(T_start, T_end, s0, lambda, alpha_beta, T_v, T_u, epsilon, h, N);
    ns(i) = length(T_history_shooting(:,1));
end
size = 20;
plot(s0s, ns, '-', 'LineWidth', 1.5, 'Color', cool(1)*0.8)
title('Shooting - Number of Iterations as a Function of $s_0$','FontSize',size, 'Interpreter','latex');
subtitle(sprintf('$N=%d$, $|$ $\\varepsilon=%g$', N, epsilon), 'Interpreter','latex')
ylabel('$n$ $[-]$','FontSize',size, 'Interpreter','latex')
xlabel('$s_0$ $[-]$','FontSize',size, 'Interpreter','latex')
grid on
grid minor
box on

% exportgraphics(fig11, 'images/shooting - n vs s0 - with negative.png','Resolution',400);

%%
% Results
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
parameters

fig12 = figure ('Position',[0 100 900 500]);

T_history = finite_difference_method(T_start, T_end, lambda, alpha_beta, T_v, T_u, epsilon, h, N);
plot(linspace(zeta_min, zeta_max, N+2), T_history{end,:}, '-', 'LineWidth', 1.5, 'Color', cool(1)*0.8)

size = 20;
title('FD - Temperature as a Function of $\zeta$','FontSize',size, 'Interpreter','latex');
subtitle(sprintf('$\\varepsilon=%g$ $|$ $N=%d$', epsilon, N), 'Interpreter','latex')
ylabel('$T$ $[-]$','FontSize',size, 'Interpreter','latex')
xlabel('$\zeta$ $[-]$','FontSize',size, 'Interpreter','latex')
grid on
grid minor
box on
% exportgraphics(fig12, 'images/FD - T vs zeta.png','Resolution',400);

%%
parameters

s0 = 0;

fig13 = figure ('Position',[150 100 900 500]);

T_history = shooting_method_with_zero(T_start, T_end, s0, lambda, alpha_beta, T_v, T_u, epsilon, h, N);
plot(linspace(zeta_min, zeta_max, N+2), T_history{end,:}, '-', 'LineWidth', 1.5, 'Color', cool(1)*0.8)

size = 20;
title('Shooting - Temperature as a Function of $\zeta$','FontSize',size, 'Interpreter','latex');
subtitle(sprintf('$\\varepsilon=%g$ $|$ $N=%d$ $|$ $s_0=%g$', epsilon, N, s0), 'Interpreter','latex')
ylabel('$T$ $[-]$','FontSize',size, 'Interpreter','latex')
xlabel('$\zeta$ $[-]$','FontSize',size, 'Interpreter','latex')
grid on
grid minor
box on
% exportgraphics(fig13, 'images/shooting - T vs zeta.png','Resolution',400);

%%
parameters

fig14 = figure ('Position',[300 100 900 500]);

T_history = finite_difference_method(T_start, T_end, lambda, alpha_beta, T_v, T_u, epsilon, h, N);
md = calc_md(T_history{end,:}, md_start, md_end, alpha_beta, lambda, h, N);
plot(linspace(zeta_min, zeta_max, N+2), md, '-', 'LineWidth', 1.5, 'Color', cool(1)*0.8)

size = 20;
title('FD - Mass Fraction as a Function of $\zeta$','FontSize',size, 'Interpreter','latex');
subtitle(sprintf('$\\varepsilon=%g$ $|$ $N=%d$', epsilon, N), 'Interpreter','latex')
ylabel('$m_d$ $[-]$','FontSize',size, 'Interpreter','latex')
xlabel('$\zeta$ $[-]$','FontSize',size, 'Interpreter','latex')
grid on
grid minor
box on
% exportgraphics(fig14, 'images/FD - md vs zeta.png','Resolution',400);

%%
parameters

s0 = 0;

fig15 = figure ('Position',[450 100 900 500]);

T_history = shooting_method_with_zero(T_start, T_end, s0, lambda, alpha_beta, T_v, T_u, epsilon, h, N);
md = calc_md(T_history{end,:}, md_start, md_end, alpha_beta, lambda, h, N);
plot(linspace(zeta_min, zeta_max, N+2), md, '-', 'LineWidth', 1.5, 'Color', cool(1)*0.8)

size = 20;
title('Shooting - Mass Fraction as a Function of $\zeta$','FontSize',size, 'Interpreter','latex');
subtitle(sprintf('$\\varepsilon=%g$ $|$ $N=%d$ $|$ $s_0=%g$', epsilon, N, s0), 'Interpreter','latex')
ylabel('$m_d$ $[-]$','FontSize',size, 'Interpreter','latex')
xlabel('$\zeta$ $[-]$','FontSize',size, 'Interpreter','latex')
grid on
grid minor
box on
% exportgraphics(fig15, 'images/shooting - md vs zeta.png','Resolution',400);

%%
parameters

s0 = 0;

fig16 = figure ('Position',[600 100 900 500]);

T_history = shooting_method(T_start, T_end, s0, lambda, alpha_beta, T_v, T_u, epsilon, h, N);
md = calc_md(T_history{end,:}, md_start, md_end, alpha_beta, lambda, h, N);
plot(linspace(zeta_min, zeta_max, N+2), md, '-', 'LineWidth', 1.5, 'Color', cool(1)*0.8)

size = 20;
title('Shooting - Mass Fraction as a Function of $\zeta$','FontSize',size, 'Interpreter','latex');
subtitle(sprintf('$\\varepsilon=%g$ $|$ $N=%d$ $|$ $s_0=%g$', epsilon, N, s0), 'Interpreter','latex')
ylabel('$m_d$ $[-]$','FontSize',size, 'Interpreter','latex')
xlabel('$\zeta$ $[-]$','FontSize',size, 'Interpreter','latex')
grid on
grid minor
box on
% exportgraphics(fig16, 'images/shooting - md vs zeta - no zero.png','Resolution',400);
% exportgraphics(fig16, 'images/shooting - md vs zeta - no zero.png','Resolution',400); exportgraphics(fig15, 'images/shooting - md vs zeta.png','Resolution',400); exportgraphics(fig14, 'images/FD - md vs zeta.png','Resolution',400); exportgraphics(fig13, 'images/shooting - T vs zeta.png','Resolution',400); exportgraphics(fig12, 'images/FD - T vs zeta.png','Resolution',400); exportgraphics(fig11, 'images/shooting - n vs s0 - with negative.png','Resolution',400); exportgraphics(fig10, 'images/shooting - n vs s0 - no negative.png','Resolution',400); exportgraphics(fig8, 'images/shooting - T vs zeta for diff epsilon.png','Resolution',400); exportgraphics(fig7, 'images/shooting - n vs epsilon.png','Resolution',400); exportgraphics(fig6, 'images/shooting - n vs N.png','Resolution',400); exportgraphics(fig5, 'images/shooting - T vs zeta for diff N.png','Resolution',400); exportgraphics(fig4, 'images/FD - T vs zeta for diff epsilon.png','Resolution',400); exportgraphics(fig3, 'images/FD - n vs epsilon.png','Resolution',400); exportgraphics(fig2, 'images/FD - n vs N.png','Resolution',400); exportgraphics(fig1, 'images/FD - T vs zeta for diff N.png','Resolution',400);








%%
% Functions
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [T_history] = finite_difference_method(T_start, T_end, lambda, alpha_beta, T_v, T_u, epsilon, h, N)
    % check that number of elements is odd
    if mod(N,2) ~= 1
        fprintf('N must be odd. It is: %d\n', N);
        T_history = 0;
        return
    end
    % calc index of zeta = 0
    index_of_zero = (N + 1) / 2 + 1; % the plus one after the fraction is because of matlab
    % set initial guess of the temperature
    T_current = T_start + (T_end - T_start) / (N+1) * [0:N+1];
    % add the current temperature into the temperature history
    T_history{1,:} = T_current;
    % the main loop of the solver
    for i = 1:1e6
        i
        % getting the updated temprature vector
        T_next = finite_difference_method_step(T_current, lambda, alpha_beta, T_v, T_u, index_of_zero, h, N);
        % add the current temperature into the temperature history
        T_history{end+1,:} = T_next;
        % check convergence
        if check_convergence_finite_difference_method(T_next, T_current, epsilon, N)
            break
        end
        % updating the temperature field
        T_current = T_next;
    end
end

function converged = check_convergence_finite_difference_method(T_next, T_current, epsilon, N)
    converged = true;
    for i = [1:N]+1
        if abs(T_next(i) - T_current(i)) > epsilon
            converged = false;
            return
        end
    end
end

function [T_next] = finite_difference_method_step(T_current, lambda, alpha_beta, T_v, T_u, index_of_zero, h, N)
    % setting boundary conditions
    T_next(1) = T_current(1);
    T_next(N+1+1) = T_current(N+1+1);
    
    % setting the temperature at zeta = 0 to be zero
    T_next(index_of_zero) = 0;

    % preforming the step
    for i = [1:N]+1
        % keeping the temperature at zeta = 0 to be zero
        if i == index_of_zero 
            continue;
        end
        T_next(i) = 0.5*(T_current(i+1) + T_current(i-1)) + 0.25 * h * lambda * exp(T_current(i)) * (T_current(i+1) - T_current(i-1)) - 0.5 * h^2 * lambda * exp(T_current(i)) * (T_v - T_u + alpha_beta);
    end
end

% ///////////////////////////////////////////////////////////////////

function [T] = shooting_method_step(T0, s0, lambda, alpha_beta, T_v, T_u, h, N)
    T(1) = T0;
    s(1) = s0;
    for i = [0:N]+1 % semi implicit Euler
        s(i+1) = lambda * exp(T(i)) * ((T_v - T_u + alpha_beta) - s(i)) * h + s(i);
        T(i+1) = s(i+1) * h + T(i);
    end
end

function [s0_next] = guess_next_s(s0_current, s0_past, F_current, F_past)
    s0_next = s0_current - F_current * (s0_current - s0_past) / (F_current - F_past);
end

function [T_history] = shooting_method(T_start, T_end, s0, lambda, alpha_beta, T_v, T_u, epsilon, h, N)
    F = [];
    s0_s = [];
    % setting initial guess
    s0_s(1) = s0;
    % the main loop of the solver
    for i = 1:1e6
        i
        % getting temperature vector according to the latest initial guess
        T = shooting_method_step(T_start, s0_s(end), lambda, alpha_beta, T_v, T_u, h, N);
        % add the current temperature into the temperature history
        if i ==1
            T_history{1,:} = T;
        else
            T_history{end+1,:} = T;
        end
        % calculating the temperature difference at the end
        F(end+1) = T(end) - T_end;
        % checking if the size of the difference is small enought
        if abs(F(end)) < epsilon
            break
        end
        % updating the initial guess of s
        if i == 1
            s0_s(end+1) = guess_next_s(s0_s(end), s0_s(end)-1, F(end), F(end)-1); % gussing initial slop of 1
        else 
            s0_s(end+1) = guess_next_s(s0_s(end), s0_s(end-1), F(end), F(end-1));
        end
    end
end

function [T_history] = shooting_method_with_zero(T_start, T_end, s0, lambda, alpha_beta, T_v, T_u, epsilon, h, N)
    % check that number of elements is odd
    if mod(N,2) ~= 1
        fprintf('N must be odd. It is: %d\n', N);
        T_history = 0;
        return
    end
    % calc index of zeta = 0
    index_of_zero = (N + 1) / 2 + 1; % the plus one after the fraction is because of matlab
    % shooting from the initial temperature to zeta = 0 and setting T = 0
    % there
    [T_history_to_zero] = shooting_method(T_start, 0, s0, lambda, alpha_beta, T_v, T_u, epsilon, h, (N+1)/2-1);
    % shooting from zeta = 0 wher T = 0 to the final temperature
    [T_history_from_zero] = shooting_method(0, T_end, s0, lambda, alpha_beta, T_v, T_u, epsilon, h, (N+1)/2-1);
    % setting the number of steps for each side
    len_to_zero = length(T_history_to_zero);
    len_from_zero = length(T_history_from_zero);
    % zipping the history of the temperature vectors according to the number of steps it
    % took
    if len_to_zero > len_from_zero
        for i = 1:len_from_zero
            T_history{i,:} = zip_two_arrays(T_history_to_zero{i,:}, T_history_from_zero{i,:}, 1);
        end
        for i = len_from_zero+1:len_to_zero
            T_history{i,:} = zip_two_arrays(T_history_to_zero{i,:}, T_history_from_zero{end,:}, 1);
        end
    elseif len_from_zero > len_to_zero
        for i = 1:len_to_zero
            T_history{i,:} = zip_two_arrays(T_history_to_zero{i,:}, T_history_from_zero{i,:}, 1);
        end
        for i = len_to_zero+1:len_from_zero
            T_history{i,:} = zip_two_arrays(T_history_to_zero{end,:}, T_history_from_zero{i,:}, 1);   
        end
    else
        for i = 1:len_to_zero
            T_history{i,:} = zip_two_arrays(T_history_to_zero{i,:}, T_history_from_zero{i,:}, 1);
        end
    end
end

function [array] = zip_two_arrays(a1, a2, num_of_overlap)
    array = [a1,a2(num_of_overlap+1:end)];
end

function md = calc_md(T, md_start, md_end, alpha_beta, lambda, h, N)
    md(1) = md_start;
    md(N+1+1) = md_end;
    for i = [1:N]+1
        md(i) = (alpha_beta * lambda * exp(T(i)))^(-1) * (T(i+1) - 2 * T(i) + T(i-1)) / h^2;
    end
end




