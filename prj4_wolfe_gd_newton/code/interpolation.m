function alpha = interpolation(x, fx, dfx, y, fy)
    a = (fy-fx-(y-x)*dfx)/(y-x)^2;
    alpha = (2*a*x-dfx)/(2*a); 
end