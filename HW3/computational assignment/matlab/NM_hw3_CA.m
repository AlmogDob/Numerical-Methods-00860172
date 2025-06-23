clc; clear; close all;

% flow_field_init = init_flow_field(ni, nj, i_AB, j_AB, i_BC, j_BC, i_CD, j_CD, i_DE, j_EF, i_GH, j_GH, i_HA, j_HA);
% figure
% contourf(x_mat, y_mat, flow_field_init, 100, "LineStyle","none")
% colormap('turbo')
% colorbar()
% axis equal

parameters
factor = 5;
% init_val = 0.1;

result = solve_Gauss_Seidel(x_min, x_max, y_min, y_max, c, mu, factor, init_val, epsilon);

fig1 = figure('Name','1','Position',[0, 250, 900, 600]);
contourf(result.x_mat, result.y_mat, result.flow_field, 200, "LineStyle","none")
colormap('turbo')
colorbar()
axis equal
title('$\phi$ Distribution Over The Hole Channel','FontSize',20,'Interpreter','latex')
xlabel('x $[in]$','FontSize',20,'Interpreter','latex')
ylabel('y $[in]$','FontSize',20,'Interpreter','latex')
% exportgraphics(fig1, 'images/phi ditribution.png','Resolution',400);


%% effect of factor =======================================================

parameters
factors = 1:1:15;
results = {};
lg = {};
for factor_index = 1:length(factors)
    factor = factors(factor_index);
    results{end+1} = solve_Gauss_Seidel(x_min, x_max, y_min, y_max, c, mu, factor, init_val, epsilon);
end

fig2 = figure('Name','2','Position',[0, 250, 900, 600]);
size = 20;
colors = cool(length(factors));
max_ni_nj_vec = [];
n_vec = [];
for factor_index = 1:length(factors)
    max_ni_nj_vec(factor_index) = max(results{factor_index}.ni, results{factor_index}.nj);
    n_vec(factor_index) = results{factor_index}.n;
end 

semilogy(max_ni_nj_vec, n_vec, 'LineStyle','-','LineWidth',2,'Color',colors(end,:))

title('Effect of Number of Elements on Number of Iterations', 'FontSize', size,'Interpreter','latex')
subtitle(sprintf('$\\varepsilon=%g$ $|$ init val: %g', results{end}.epsilon, results{end}.init_val), 'FontSize',size-4,'Interpreter','latex')
ylabel('Number of Iterations', 'FontSize', size,'Interpreter','latex')
xlabel('$\max{\left(ni,nj\right)}$', 'FontSize', size,'Interpreter','latex')
grid on 
grid minor
box on
% exportgraphics(fig2, 'images/ni nj - n.png','Resolution',400);

fig3 = figure('Name','3','Position',[200, 250, 900, 600]);
size = 20;
colors = cool(length(factors));
max_vec = [];
for factor_index = 1:length(factors)-1
    max_vec(factor_index) = calc_diff(results{factor_index}, results{factor_index+1});
end 

semilogy(max_ni_nj_vec(1:end-1), max_vec, 'LineStyle','-','LineWidth',2,'Color',colors(end,:))

title('Effect of Number of Elements on Maximum Value', 'FontSize', size,'Interpreter','latex')
subtitle(sprintf('$\\varepsilon=%g$ $|$ init val: %g', results{end}.epsilon, results{end}.init_val), 'FontSize',size-4,'Interpreter','latex')
ylabel('$\left|\phi_{max}^{n+1}-\phi_{max}^{n}\right|$', 'FontSize', size,'Interpreter','latex')
xlabel('$\max{\left(ni,nj\right)}$', 'FontSize', size,'Interpreter','latex')
grid on 
grid minor
box on
% exportgraphics(fig3, 'images/ni nj - max diff.png','Resolution',400);
% exportgraphics(fig2, 'images/ni nj - n.png','Resolution',400); exportgraphics(fig3, 'images/ni nj - max diff.png','Resolution',400);


%% effect of epsilon ======================================================

parameters
epsilons = logspace(-1, -17, 17);
% epsilons = logspace(0, -5, 6);
results = {};
lg = {};
for epsilons_index = 1:length(epsilons)
    epsilon = epsilons(epsilons_index);
    results{end+1} = solve_Gauss_Seidel(x_min, x_max, y_min, y_max, c, mu, factor, init_val, epsilon);
end

fig4 = figure('Name','4','Position',[400, 250, 900, 600]);
size = 20;
colors = cool(length(epsilons));
max_ni_nj_vec = [];
n_vec = [];
for epsilons_index = 1:length(epsilons)
    n_vec(epsilons_index) = results{epsilons_index}.n;
end 

semilogx(epsilons, n_vec, 'LineStyle','-','LineWidth',2,'Color',colors(end,:))

title('Effect of $\phi$ on Number of Iterations', 'FontSize', size,'Interpreter','latex')
subtitle(sprintf('$\\max{\\left(ni,nj\\right)}=%d$ $|$ init val: %g', max(results{end}.ni, results{end}.nj), results{end}.init_val), 'FontSize',size-4,'Interpreter','latex')
ylabel('Number of Iterations', 'FontSize', size,'Interpreter','latex')
xlabel('$\varepsilon$', 'FontSize', size,'Interpreter','latex')
grid on 
grid minor
box on
% exportgraphics(fig4, 'images/epsilon - n.png','Resolution',400);

fig5 = figure('Name','5','Position',[600, 250, 900, 600]);
size = 20;
colors = cool(length(epsilons));
max_vec = [];
for epsilons_index = 1:length(epsilons)-1
    max_vec(epsilons_index) = calc_diff(results{epsilons_index}, results{epsilons_index+1});
end 

loglog(epsilons(1:end-1), max_vec, 'LineStyle','-','LineWidth',2,'Color',colors(end,:))

title('Effect of $\phi$ on Maximum Value', 'FontSize', size,'Interpreter','latex')
subtitle(sprintf('$\\max{\\left(ni,nj\\right)}=%d$ $|$ init val: %g', max(results{end}.ni, results{end}.nj), results{end}.init_val), 'FontSize',size-4,'Interpreter','latex')
ylabel('$\left|\phi_{max}^{n+1}-\phi_{max}^{n}\right|$', 'FontSize', size,'Interpreter','latex')
xlabel('$\varepsilon$', 'FontSize', size,'Interpreter','latex')
grid on 
grid minor
box on
% exportgraphics(fig5, 'images/epsilon - max diff.png','Resolution',400);
% exportgraphics(fig5, 'images/epsilon - max diff.png','Resolution',400); exportgraphics(fig4, 'images/epsilon - n.png','Resolution',400);


%% effect of init val =====================================================

parameters
init_vals = logspace(-2,4,50);
results = {};
lg = {};
for init_vals_index = 1:length(init_vals)
    init_val = init_vals(init_vals_index);
    results{end+1} = solve_Gauss_Seidel(x_min, x_max, y_min, y_max, c, mu, factor, init_val, epsilon);
end

fig6 = figure('Name','6','Position',[0, 150, 900, 600]);
size = 20;
colors = cool(length(init_vals));
max_ni_nj_vec = [];
n_vec = [];
for init_vals_index = 1:length(init_vals)
    n_vec(init_vals_index) = results{init_vals_index}.n;
end 

loglog(init_vals, n_vec, 'LineStyle','-','LineWidth',2,'Color',colors(end,:))

title('Effect of Initial values on Number of Iterations', 'FontSize', size,'Interpreter','latex')
subtitle(sprintf('$\\max{\\left(ni,nj\\right)}=%d$ $|$ $\\varepsilon=%g$', max(results{end}.ni, results{end}.nj), results{end}.epsilon), 'FontSize',size-4,'Interpreter','latex')
ylabel('Number of Iterations', 'FontSize', size,'Interpreter','latex')
xlabel('Initial values', 'FontSize', size,'Interpreter','latex')
grid on 
grid minor
box on
% exportgraphics(fig6, 'images/init val - n.png','Resolution',400);

fig7 = figure('Name','7','Position',[200, 150, 900, 600]);
size = 20;
colors = cool(length(init_vals));
max_vec = [];
for init_vals_index = 1:length(init_vals)-1
    max_vec(init_vals_index) = calc_diff(results{init_vals_index}, results{init_vals_index+1});
end 

loglog(init_vals(1:end-1), max_vec, 'LineStyle','-','LineWidth',2,'Color',colors(end,:))

title('Effect of Initial values on Maximum Value', 'FontSize', size,'Interpreter','latex')
subtitle(sprintf('$\\max{\\left(ni,nj\\right)}=%d$ $|$ $\\varepsilon=%g$', max(results{end}.ni, results{end}.nj), results{end}.epsilon), 'FontSize',size-4,'Interpreter','latex')
ylabel('$\left|\phi_{max}^{n+1}-\phi_{max}^{n}\right|$', 'FontSize', size,'Interpreter','latex')
xlabel('Initial values', 'FontSize', size,'Interpreter','latex')
grid on 
grid minor
box on
% exportgraphics(fig7, 'images/init val - max diff.png','Resolution',400);
% exportgraphics(fig7, 'images/init val - max diff.png','Resolution',400); exportgraphics(fig6, 'images/init val - n.png','Resolution',400);





%% Functions ==============================================================

function flow_field = set_BC(flow_field, i_AB, j_AB, i_BC, j_BC, i_CD, j_CD, i_DE, j_EF, i_GH, j_GH, i_HA, j_HA)
    for i = i_AB+1
        for j = j_AB+1
            flow_field(j, i) = 0;
        end
    end
    for i = i_BC+1
        for j = j_BC+1
            flow_field(j, i) = 0;
        end
    end
    for i = i_CD+1
        for j = j_CD+1
            flow_field(j, i) = 0;
        end
    end
    for i = i_DE+1
        for j = j_EF+1
            flow_field(j, i) = 0;
        end
    end
    for i = i_GH+1
        for j = j_GH+1
            flow_field(j, i) = 0;
        end
    end
    for i = i_HA+1
        for j = j_HA+1
            flow_field(j, i) = 0;
        end
    end
end

function flow_field = init_flow_field(ni, nj, init_val, i_AB, j_AB, i_BC, j_BC, i_CD, j_CD, i_DE, j_EF, i_GH, j_GH, i_HA, j_HA)
    flow_field = ones(nj+1,ni+1) * init_val;
    flow_field = set_BC(flow_field, i_AB, j_AB, i_BC, j_BC, i_CD, j_CD, i_DE, j_EF, i_GH, j_GH, i_HA, j_HA);
end

function converged = check_convergence(flow_field_next, flow_field_current, epsilon,j_AB, i_BC)
    converged = true;
    for i = i_BC+1
        for j = j_AB+1
            if abs(flow_field_next(j, i) - flow_field_current(j, i)) > epsilon
                converged = false;
                return
            end
        end
    end
end

function result = solve_Gauss_Seidel(x_min, x_max, y_min, y_max, c, mu, factor, init_val, epsilon)
    calc_extra_param
    flow_field_current = init_flow_field(ni, nj, init_val, i_AB, j_AB, i_BC, j_BC, i_CD, j_CD, i_DE, j_EF, i_GH, j_GH, i_HA, j_HA);
    flow_field_next = flow_field_current;
    for n = 1:1e6
        if ~mod(n,100)
            fprintf('factor: %d | epsilon: %g | init val: %g | n: %d\n', factor, epsilon, init_val, n)
        end
        for i = i_BC(2:end-1)+1
            for j = j_AB(2:end-1)+1
                % flow_field_next(j, i) = 0.5 * ((flow_field_next(j, i-1) + flow_field_current(j, i+1)) * delta_y^2 / (delta_y^2+delta_x^2) + (flow_field_next(j-1, i) + flow_field_current(j+1, i)) * delta_x^2 / (delta_y^2+delta_x^2) + delta_x^2*delta_y^2 / (delta_y^2+delta_x^2) * c / mu);
                flow_field_next(j, i) = 0.5 * ((flow_field_current(j, i-1) + flow_field_current(j, i+1)) * delta_y^2 / (delta_y^2+delta_x^2) + (flow_field_current(j-1, i) + flow_field_current(j+1, i)) * delta_x^2 / (delta_y^2+delta_x^2) + delta_x^2*delta_y^2 / (delta_y^2+delta_x^2) * c / mu);
            end
        end
        flow_field_next = set_BC(flow_field_next, i_AB, j_AB, i_BC, j_BC, i_CD, j_CD, i_DE, j_EF, i_GH, j_GH, i_HA, j_HA);
        if check_convergence(flow_field_next, flow_field_current, epsilon,j_AB, i_BC)
            break
        end
        flow_field_current = flow_field_next;
    end
    result.n          = n;
    result.x_min      = x_min;
    result.x_max      = x_max;
    result.y_min      = y_min;
    result.y_max      = y_max;
    result.c          = c;
    result.mu         = mu;
    result.factor     = factor;
    result.epsilon    = epsilon;
    result.x_mat      = x_mat;
    result.y_mat      = y_mat;
    result.delta_x    = delta_x;
    result.delta_y    = delta_y;
    result.ni         = ni;
    result.nj         = nj;
    result.i_BC       = i_BC;
    result.j_AB       = j_AB;
    result.init_val   = init_val;
    result.flow_field = flow_field_next;
end

% function rms = RMS(result1, result2)
% % This functions calcutlates the RMS value for 1 compared to 2
%     flow_field2_at_coord_of_1 = [];
%     for i = result1.i_BC+1
%         for j = result1.j_AB+1
%             flow_field2_at_coord_of_1(j, i) = interp2(result2.x_mat, result2.y_mat, result2.flow_field, result1.x_min + i*result1.delta_x, result1.y_min + j*result1.delta_y);
%         end
%     end
%     flow_field2_at_coord_of_1(isnan(flow_field2_at_coord_of_1)) = 0;
%     flow_field1 = result1.flow_field;
%     rms = sqrt(sum(sum((flow_field1-flow_field2_at_coord_of_1).^2)) / (length(flow_field1(:,1)) * length(flow_field1(1,:))));
% end

function diff = calc_diff(result1, result2)
    max1 = max(max(result1.flow_field));
    max2 = max(max(result2.flow_field));
    diff = abs(max1-max2);
end
