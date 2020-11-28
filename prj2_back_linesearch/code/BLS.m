function data = BLS(f, step, rho, c, x0, maxiter, tol, mode)
    % data = BLS(f, step, rho, c, x0, maxiter, tol)
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
    % mode: 0-1 variable, 0-steepest descend, 1-newton method.
    fc = feval(f, x0(1), x0(2));
    gc = Diff(x0(1), x0(2));
    iter = 0;
    err = 1.0;
    x = x0;
    while iter < maxiter && abs(err) > tol
        iter = iter + 1;        
        alpha = step;
        gc = Diff(x(1), x(2));
        B = Hessian(x(1), x(2));
        if mode == 1 % newton or steepest descend
            pk = -inv(B)*gc; 
        else
            pk = -gc/norm(gc);
        end
        xn = x + alpha * pk';
        while feval(f, xn(1), xn(2)) > fc + c * alpha * gc'*pk
            alpha = rho * alpha;
            xn = x + alpha * pk';
        end
        % Update information
        fc = feval(f, xn(1), xn(2));
        x = xn;
        err = fc;
        % Record and Return
        err_rec(iter) = err;
        alpha_rec(iter) = alpha;
        x_rec(iter, :) = x(:);
    end
    data = [(1: iter)', alpha_rec', x_rec, err_rec']; % Return all the information of iterations
end

%% ------------------------------------------------------
function gradient = Diff(x1, x2)
    % gradient = Diff(x1, x2)
    % Get the Gradient vector of Rosenbrock funtion.
    % For a more common use, we can use symbolic calculation to give out
    % the corresponding vector that user requires on every iteration.
    % Deprecated: syms x1 x2, input f, xc
    % gradient = double([vpa(subs(diff(f(x1, x2), x1), [x1, x2], xc));
    %           vpa(subs(diff(f(x1, x2), x2), [x1, x2], xc))]);
    % This could calculate all kinds of two-variable functions.
    % x = (x1, x2) x1, x2: 1-d respectively
    % gradient: 2-d vector

    gradient = [(2*x1 - 400*x1*(- x1^2 + x2) - 2);  
                (- 200*x1^2 + 200 * x2)];
end

%% ----------------------------------------------
function Bk = Hessian(x1, x2)
    % Bk = Hessian(x1, x2)
    % Get the Hessian Matrix of Rosenbrock funtion.
    % For a more common use, we can use symbolic calculation to give out
    % the corresponding vector that user requires on every iteration.
    % Deprecated: syms x1 x2, input f, xc
    % gradient = double[subs(diff(diff(f(x1, x2), x1), x1), [x1, x2], xc), 
    %                   subs(diff(diff(f(x1, x2), x1), x2), [x1, x2], xc); 
    %                   subs(diff(diff(f(x1, x2), x2), x1), [x1, x2], xc), 
    %                   subs(diff(diff(f(x1, x2), x2), x2), [x1, x2], xc),
    % This could calculate all kinds of two-variable functions.
    % x = (x1, x2) x1, x2: 1-d respectively
    % gradient: 2-d vector
    
    Bk = [1200*x1^2 - 400*x2 + 2, -400*x1; 
          -400*x1,  200];
end