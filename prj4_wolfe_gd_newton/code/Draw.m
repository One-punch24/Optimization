% This script is to draw the effect of Wolfe Linesearch.
clc; clear all;
data = Wolfe(1, 1e-4, 20*rand(1,5)-10, 20000, 1e-4, 1)

plot(data(:, end),'linewidth', 2);
xlabel('Iteration');
ylabel('Error');
title('Newton');
%%
figure;
data = Wolfe(1, 1e-4, rand(1,10), 20000, 1e-4, 0)

plot(data(:, end),'linewidth', 2);
xlabel('Iteration');
ylabel('Error');
title('Gradient Descendent');