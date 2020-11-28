function [L, tau] = modifiedCholesky(A) 
    tau = 0;
    [m,n] = size(A) ; 
    if (m ~= n)
        error("Matrix must be square. Size %dx%d",m,n);
    end
    mindiag = min(diag(A));
    beta = 1e-3;
    if (mindiag > 0)
        tau = 0;
        R = A;
    else
        tau = beta - mindiag;
        A_tmp = A + tau * eye(n,n);
    end
    L = zeros(n, n);
    while 1
        A_tmp = A + tau * eye(n);
        flag = 0;
        for i = 1: n
            if i == 1
                if A_tmp(i, i) > 0
                    L(i, i) = sqrt(A_tmp(i, i));
                    L(:, i) = A_tmp(:, i)./L(i, i);
                else
                    flag = 1;
                    break;
                end
            else
                if A_tmp(i, i)-sum(L(i, 1:i-1).^2) > 0
                    L(i, i) = sqrt(A_tmp(i, i)-sum(L(i, 1:i-1).^2));
                    if i < n
                        s = zeros(n-i, 1);
                        for ind = 1: i-1
                            s = s + L(i+1: n, ind).*L(i, ind);
                        end
                        L(i+1:n, i) = (A_tmp(i+1:n, i)-s)./L(i, i);
                    end
                else
                    flag = 1;
                    break;
                end
            end     
        end
        if flag == 0
            break;
        else
            tau = max([2 * tau, beta]);
        end
    end
end