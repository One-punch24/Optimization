function [xn, delta, tau] = trustregion(Delta, delta, x, fk, gk, Bk, t, y)
% [xn, delta] = trustregion(Delta, delta, x, fk, gk, Bk, t, y)
    gk
    Bk
    eta = 1e-4;
%     delta = Delta;
    mk = @(p) fk+gk'*p+0.5*p'*Bk*p;
    pn = -Bk\gk;
    pu = -(gk'*gk)/(gk'*Bk*gk)*gk;
    tau = dogleg(pn, pu, delta);
    if tau < 1
        pk = tau * pu;
    else
        pk = pu + (tau-1)*(pn-pu);
    end
    rho = (fk - f(x+pk, t, y))/(mk(zeros(5, 1))-mk(pk));
    if rho < 1/4
        delta = delta/4.0;
    else
        if rho > 3/4 && norm(pk) == delta
            delta = min(2*delta, Delta);
        else
            delta = delta;
        end
    end
    if rho > eta
        xn = x + pk;
    else
        xn = x;
    end
     
end

function tau = dogleg(pn, pu, delta)
    if norm(pn) < delta
        tau = 2;
    else
        if norm(pu) >= delta % tau <= 1
            tau = delta / norm(pu);
        else
            a = norm(pn-pu)^2;
            b = 2*(2*pu-pn)'*(pn-pu);
            c = norm(2*pu-pn)^2-delta^2;
            tau = (-b+sqrt(b^2-4*a*c))/(2*a);
            if tau < 1 || tau > 2
                tau = (-b-sqrt(b^2-4*a*c))/(2*a);
            end
        end
    end    
%     tau
end
function fx = f(x, t, y)
    C1 = exp(-x(1));
    C2 = exp(-x(2));
    C3 = exp(-x(3));
    C4 = x(4);
    C5 = x(5);
    res = C5*C1.^(C2.^((1-C3.^t).^C4))-y;
    fx = 0.5 * norm(res, 2)^2;
end