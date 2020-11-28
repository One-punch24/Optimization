function B = LDLT(A,delta)
    
    if ~ishermitian(A)
        error('Pleae input a symmetric matrix.')
    end
    if nargin < 2
        delta = 1e-3;
    end

    n = max(size(A));
    A0 = A;
    k = 1;
    DMC = eye(n);
    D = eye(n);
    L = eye(n);
    pp = 1:n;
    normA = norm(A(:),inf);
    rho = normA;

    alpha = 0.64; %(1 + sqrt(17))/8;

    while k < n
          [lambda, r] = max(abs(A(k+1:n,k)));
          if lambda > 0
              if abs(A(k,k)) >= alpha*lambda
                  s = 1;
              else
                  j = k;
                  pivot = 0;
                  lambda_j = lambda;
                  while ~pivot
                        [temp_, r] = max(abs(A(k:n,j)));
                        r = r(1) + k - 1; % If more than one, take the smallest index.
                        temp = A(k:n,r); temp(r-k+1) = 0;
                        lambda_r = max(abs(temp));
                        if alpha*lambda_r <= abs(A(r,r))
                           pivot = 1;
                           s = 1;                       
                           A([k, r], :) = A([r, k], :);
                           L([k, r], :) = L([r, k], :);
                           A(:, [k, r]) = A(:,[r, k]);
                           L(:, [k, r]) = L(:,[r, k]);
                           pp([k, r]) = pp([r, k]);
                        elseif lambda_j == lambda_r
                           pivot = 1;
                           s = 2;
                           % Swap the cols and rows.
                           A([k, j], :) = A([j, k], :);
                           L([k, j], :) = L([j, k], :);
                           A(:, [k, j]) = A(:, [j, k]);
                           L(:, [k, j]) = L(:, [j, k]);
                           pp([k, j]) = pp([j, k]);

                           A([k+1, r], :) = A([r, k+1], :);
                           L([k+1, r], :) = L([r, k+1], :);
                           A(:, [k+1, r]) = A(:, [r, k+1]);
                           L(:, [k+1, r]) = L(:, [r, k+1]);
                           pp([k+1, r]) = pp([r, k+1]);
                        else
                           j = r;
                           lambda_j = lambda_r;               
                        end
                  end
              end

              if s == 1

                 D(k,k) = A(k,k);
                 A(k+1:n,k) = A(k+1:n,k)/A(k,k);
                 L(k+1:n,k) = A(k+1:n,k);
                 i = k+1:n;
                 A(i,i) = A(i,i) - A(i,k) * A(k,i);
                 A(i,i) = 0.5 * (A(i,i) + A(i,i)');

              elseif s == 2

                 E = A(k:k+1,k:k+1);
                 D(k:k+1,k:k+1) = E;
                 index = k+2:n;
                 C = A(index ,k:k+1);
                 temp = C/E;
                 L(index ,k:k+1) = temp;
                 A(index, index) = A(index, index) - temp*C';
                 A(index, index) = 0.5 * (A(index, index) + A(index, index)');

              end
          else  % Nothing to do.

             s = 1;
             D(k,k) = A(k,k);

          end
          % Modified Cholesky
          if s == 1
             DMC(k, k) = max(D(k, k), delta);
          elseif s == 2 
             E = D(k:k+1,k:k+1);
             [U,T] = eig(E); % Orthogonal Decomposition
             for ii = 1:2
                 T(ii, ii) = max(T(ii, ii), delta);
             end
             temp = U*T*U';
             DMC(k:k+1,k:k+1) = (temp + temp')/2;  % Ensure symmetric.

          end
          k = k + s;
          if k == n
             D(n,n) = A(n,n);
             DMC(k, k) = max(D(k, k), delta);
             break
          end
    end
    if min(min(D)) < 0
        B=L*DMC*L';
    else
        B = A0;
    end
end
