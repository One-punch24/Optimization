%% Load Dataset
clear all;
load dataset1.txt
global t y
t = dataset1(:, 1); % Variable
y = dataset1(:, 2); 
%% Iterations
X = [10; 10; 0.1; 10; 150];   % Param Initialization.
Delta = 10;
delta = Delta*0.8;
rec = [];
rec_tau=[];
for iter = 1: 30
    [fx, J] = jacobi(t, y, X);
    H = J'*J;
    H = LDLT(H, 1e-3);
    B = -J'* fx;
    rec = [rec; norm(fx)];    % record the norm of error vector.
    [Xn, delta, tau] = trustregion(Delta, delta, X, 0.5*norm(fx,2)^2, -B, H, t, y);
    rec_tau = [rec_tau; tau]; % record the tau value
    if norm(fx) < 1
        break
    end
    X = Xn;
end
visualization(X);
figure; plot(rec, 'linewidth', 2); xlabel('Iterations'); ylabel('Error');

%% Supportive Function
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