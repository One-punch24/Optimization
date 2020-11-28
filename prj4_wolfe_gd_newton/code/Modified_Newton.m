function data = Modified_Newton(step, c, x0, maxiter, tol, mode)
    % data = BLS(step, rho, c, x0, maxiter, tol, mode)
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
    % mode: 0-no add, 1-chol. modified newton method.

    f = @func;
    fc = feval(f, x0);
    gc = dfunc(x0);
    iter = 0;
    err = 1.0;
    x = x0;
    alpha = step;
    c2 = 0.9;

    while iter < maxiter && abs(err) > tol
        iter = iter + 1;        
             
        alpha = step;
        a = 0; b = Inf; 
        gc = dfunc(x);
        B = d2func(x);
        if mode == 1
            B = LDLT(B, 1e-3);
%           [L, tau]  modifiedCholesky(B);
%            B = L*L';
        end
        pk = -B\gc;
        xn = x + alpha * pk';
        while 1
            if func(xn) > fc + c * alpha * gc'*pk
                b = alpha;
                alpha = (alpha+a)/2;
                xn = x + alpha * pk';
                continue;
            end
            if dfunc(xn)'*pk < c2*dfunc(x)'*pk
                a = alpha;
                alpha = min([2*alpha, (b+alpha)/2]);
                xn = x + alpha * pk';
                if abs(a-b) < 1e-6
                    break
                end
                continue;
            end
            break;
        end
      % alpha
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
