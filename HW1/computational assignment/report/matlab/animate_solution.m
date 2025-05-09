clc; clear; close all;

parameters

T_history_FD = finite_difference_method(T_start, T_end, lambda, alpha_beta, T_v, T_u, epsilon, h, N);
T_history_S = shooting_method(T_start, T_end, 0, lambda, alpha_beta, T_v, T_u, epsilon, h, N);
% [T_history_S] = shooting_method_with_zero(T_start, T_end, 0, lambda, alpha_beta, T_v, T_u, epsilon, h, N);

fig1 = figure ('Name', '1', 'Position',[0 150 900 500]);

for i = 1:1:length(T_history_S(:,1))
    i
    clf(fig1);
    hold all
    plot(linspace(zeta_min, zeta_max, N+2), T_history_FD{end,:}, '--')
    plot(linspace(zeta_min, zeta_max, N+2), T_history_S{i,:})

    drawnow
    input('press');
end









% Functions
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [T_history] = finite_difference_method(T_start, T_end, lambda, alpha_beta, T_v, T_u, epsilon, h, N)
    % check that number of elements is odd
    if mod(N,2) ~= 1
        fprintf("N must be odd. It is: %d\n", N);
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
        fprintf("N must be odd. It is: %d\n", N);
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






