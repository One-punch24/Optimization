function [s, record] = TSPSA(cityXY, T0, a, MAXIT, N)
    k = 0;
    t = T0;
    x = cityXY;
    n = length(x);
    xs = x;
    ds = distance(xs);
    record = zeros(1, N);
    for i = 1:N
        while (k < MAXIT)
            dx = distance(x);
            if (dx < ds)
                xs = x;
                %ds = distance(xs);
                ds = dx;
            end
            y = swapcities(x);
            dy = distance(y);
            h = min(1, exp(-(dy - dx)/t));
            U = rand();
            if (U < h)
                x = y;
            end
            k = k + 1;
        end
%         x = xs;
        distance(x)
        t = t * a;
%         t = t * (1 - 1 / (50 + log(10 + log(i))))
%         t = t * (1 - 1 / (50 + log(i)))
        k = 0;
        record(i) = distance(x);
    end
    s = xs;
end


