% Initial Value Dependency
count = 0;
testnum = 10000;
converge = 50;

for i = 2: 10
    num = i;
    count = 0;
    x0 = [converge*rand(1000, 1), converge*rand(1000, num-1)-converge/2];
    for iter = 1: testnum
        data = Modified_Newton(1, 1e-4, x0(iter, :), 2000, 1e-4, 1);
        length(data);
        iter
        if length(data) == 2000
            count = count + 1
            iter
        end
    end
    rate(i) = count / testnum
    pause
end

%% Convergence Rate
sumup = 0;
count = 0;
testnum = 10000;
for vari = 2: 10
    count = 0;
    sumup = 0;
%     x0 = rand(testnum, vari)+0.5 % 0.5~ 1.0
    x0 = rand(testnum, vari)*10-5; % 0.5~ 1.0
    for i = 1: testnum
        data = Modified_Newton(1, 1e-4, x0(i, :), 2000, 1e-4, 1);
        err0 = data(:, end);
        if length(data) < 2000
            index = min(find(err0 < 1));
            err = err0(max(end-20, index): end);
            est_rate = abs(mean(log(err(2:end)./err(1:end-1))));
            sumup = sumup+est_rate;
            count = count + 1
        end
    end
    i
    result = sumup/count
    pause
end