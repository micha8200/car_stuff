function [A_opt, c_opt] = fit_piecewise_model(x_data, y_data)
    % Objective function: computes residuals
    function err = residuals(params)
        A = params(1);
        c = params(2);
        y_pred = arrayfun(@(x) model_func(x, A, c), x_data);
        err = y_pred - y_data;
    end

    % Piecewise model function
    function y = model_func(x, A, c)
        if x > c
            y = A / x;
        else
            y = A / c;
        end
    end

    % Initial guess for [A, c]
    A0 = 1;
    c0 = median(x_data);
    initial_guess = [A0, c0];

    % Lower bounds to avoid division by zero
    lb = [1e-8, 1e-8];
    ub = [Inf, Inf];

    % Perform least squares optimization
    options = optimoptions('lsqnonlin', 'Display', 'off');
    [opt_params, ~] = lsqnonlin(@residuals, initial_guess, lb, ub, options);

    % Extract optimal parameters
    A_opt = opt_params(1);
    c_opt = opt_params(2);
end