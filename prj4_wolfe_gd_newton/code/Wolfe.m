function data = Wolfe(step, c, x0, maxiter, tol, mode)
    % data = LS_inter(step, rho, c, x0, maxiter, tol, mode)
    % Backtracking Line Search algorithm to adjust alpha that satisfies
    % Wolfe Conditions.
    % f: two-variable function handle
    % step: guessed alpha_0, usually 1
    % rho: the scaling factor of alpha, 0 < rho < 1
    % c: learning rate, usually 1e-4
    % x0: starting point
    % maxiter: the maximum time of iteration, usually 20
    % tol: tolerance of error, usually 1e-7
    % data: a (iter, 6) matrix consisting iteration, alpha, x, error information.
    % mode-0: GD; -1: Newton.
    fc = func(x0);
    gc = dfunc(x0);
    iter = 0;
    err = 1.0;
    x = x0;
    alpha = step;
    if mode == 1
        c2 = 0.9;
    else
        c2 = 0.1;
    end

    while iter < maxiter && abs(err) > tol
        iter = iter + 1;           
        if mode == 1
            alpha = step;
        else
            alpha = alpha*1.01;
        end
        gc = dfunc(x);
        B = d2func(x);
        if mode == 1
            pk = -B\gc;
        else
            pk = -gc/norm(gc);
        end
        xn = x + alpha * pk';
        
        % If not satisfying wolfe condition, do interpolation.
        if func(xn) > fc + c * alpha * gc'*pk || dfunc(xn)'*pk < c2*dfunc(x)'*pk
            alpha = interpolation(0, fc, gc'*pk, alpha, func(x+alpha*pk'));        
        end
        xn = x + alpha * pk';
        % Update information
        fc = func(xn);
        x = xn;
        err = norm(x-ones(size(x)));
        % Record and Return
        err_rec(iter) = err;
        alpha_rec(iter) = alpha;
        x_rec(iter, :) = x(:);
    end
    data = [(1: iter)', alpha_rec', x_rec, err_rec']; % Return all the information of iterations
end
