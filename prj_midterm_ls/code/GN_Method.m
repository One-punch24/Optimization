%% Load Dataset
clear all;
load dataset1.txt
global t y
t = dataset1(:, 1); % Variable
y = dataset1(:, 2); 

%% Iterations
X = [10; 10; 0.02; 1; 100];   % Param
rec = [];

for iter = 1: 150
    [fx, J] = jacobi(t, y, X);
    H = J'*J;
    H = LDLT(H, 1e-3);
    B = -J'* fx;
    dX = H\B;
    iter    
    rec = [rec; norm(fx)];
    alpha = linesearch(X, dX);
    Xn = X + alpha * dX;
    if norm(fx) < 1
        break
    end
    X = Xn;
end
visualization(X);
figure; plot(rec, 'linewidth', 2); xlabel('Iterations'); ylabel('Error');

%% Supportive Function
function alpha = linesearch(x, pk)
    global y t
    alpha = 1e-5; 
    c1 = 1e-4; c2 = 0.9;
    a = 0; b = Inf;
    [fx, J] = jacobi(t,y, x);
    fc = 0.5 * norm(fx, 2)^2;
    gc = J'*fx;
    xn = x + alpha * pk;
    while 1
        [fx, J] = jacobi(t,y,xn);
        gx = J'*fx;
        if fx > fc + c1 * alpha * gc'*pk
            b = alpha;
            alpha = (alpha+a)/2;
            xn = x + alpha * pk;
            continue;
        end
        if gx'*pk < c2*gc'*pk
            a = alpha;
            alpha = min([2*alpha, (b+alpha)/2]);
            xn = x + alpha * pk;
            if abs(a-b) < 1e-8
                break
            end
            continue;
        end
        break;
    end
%     alpha = min(alpha, 0.1);
end

function visualization(x)
    global y t
    C1 = exp(-x(1));
    C2 = exp(-x(2));
    C3 = exp(-x(3));
    C4 = x(4);
    C5 = x(5);
    yy = C5.*C1.^(C2.^((1-C3.^t).^C4));
    figure; plot(t, yy,  'linewidth', 2); hold on;
    plot(t, y, '*-', 'linewidth', 2); legend('Fit', 'Original');
    xlabel('PO2'); ylabel('SO2');
    
    figure; plot(t, yy-y, 'linewidth', 2);
    xlabel('PO2'); ylabel('Error of SO2');
%     plot(t, yy-y)
end