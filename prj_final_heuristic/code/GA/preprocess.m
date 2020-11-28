function D = preprocess(cityXY)
    N = length(cityXY);
    D = zeros(N, N);
    for i = 1: N
        D(i, :) = sqrt((cityXY(1, i)-cityXY(1,:)).^2+(cityXY(2, i)-cityXY(2,:)).^2);
    end
end