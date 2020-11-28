function visualize(f, data, textmode)
    % visualize(f, data) 
    % To visualize the contour map and error vs. iterations.
    % f: two-variable function handle, @rosenbrock
    % data: n*5-d, with [iterations, alpha, x1, x2, error]
    % textmode: 1-text on; 0-text off.
    
    iter = data(:, 1);
    alpha = data(:, 2);
    x = data(:, 3:4);
    error = data(:, 5);
    figure;  subplot(121); contour_r(f);
    % Adjust the visual field according to the iterative points.
    minx = min(x(:,1)); maxx = max(x(:,1)); miny = min(x(:,2)); maxy = max(x(:,2));
    scatter(x(:,1), x(:, 2), 30, 'r', 'filled');
    plot(x(:,1), x(:,2), 'r', 'linewidth', 1.5);
    if textmode == 1
        for i = 1: length(x)
            text(x(i, 1), x(i, 2), num2str(i), 'Fontsize', 20, 'Fontname', 'Consolas');
        end    
    end
    % Visualization Modules
    axis([minx-0.2, maxx+0.2, miny-0.2, maxy+0.2]);
    xlabel('x1'); ylabel('x2'); title('Contour Map Visualization');
    subplot(122);
    plot(iter, error,'linewidth',2);
    xlabel('Iterations'); ylabel('Error'); title('Error Evaluation');
end
%% ------------------------------------------------------------
function contour_r(f)
    % contour_r(r) is to draw the contour map with logrithmic intervals for
    % given function f within the certain range around (1, 1).
    % Input: f: function handle
    % Output: Contour Map, no return.
    
    x = -3.5:0.01:3.5; 
    y = -3.5:0.01:3.5;
    [X, Y] = meshgrid(x, y);
    Z = f(X, Y);
    v = logspace(-10, 2, 100);
    contour(X, Y, Z, v,  'LineWidth', 1.2); hold on;
    scatter(1, 1, 200, 'r', '*');
    axis equal;
end