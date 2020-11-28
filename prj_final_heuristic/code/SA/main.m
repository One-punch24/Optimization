clear all

%cityXY = rand(2, 31);
%x = TSPSA(cityXY, 100, 0.5, 4000, 20);

% cityXY = rand(2, 51);
%x = TSPSA(cityXY, 100, 0.5, 4000, 20);

%cityXY = rand(2, 71);
%x = TSPSA(cityXY, 100, 0.5, 4000, 20);

% cityXY = rand(2, 101);
% x = TSPSA(cityXY, 100, 0.5, 4000, 20);

% cityXY = rand(2, 101);
% x = TSPSA(cityXY, 100, 0.5, 4000, 40);
load x50.mat
% load x50.mat
xs = cityXY;
ds = 100;
iters = [500, 1000, 2000];
num_t = 200;
N = 1000;
M = length(iters);
record = [];

dl = zeros(M,N);
figure
for j = 1:1
    for i = 1:N
        i
        [x, record] = TSPSA(cityXY, 100, 0.95, iters(j), num_t);
        distance(x)
        dl(j, i) = distance(x);
        if (dl(j, i) < ds)
            xs = x;
            ds = dl(j, i);
        end
        pause
        if mod(i, 20) == 0
            hist(dl(j,1:i));
            pause(0.001)
        end
    end
    
    j
end
%%
hold off;
plotcities(xs);
axis equal;
%%
figure;
subplot(311); hist(dl(1,:)); title('500 iters');
subplot(312); hist(dl(2,:)); ylabel('Frequency');title('1000 iters');
subplot(313); hist(dl(3,:)); title('2000 iters');
% save dl.mat dl