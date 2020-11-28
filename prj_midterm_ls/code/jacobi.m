function [res, J] = jacobi(data, y, param)
%     df/dx
%     fval
    m = length(data);
    n = length(param);
    J = zeros(m, n);
%     H = zeros(n, n);
%     tmp = zeros(n, n);
    res = zeros(m, 1);
    x1 = param(1);
    x2 = param(2);
    x3 = param(3);
    x4 = param(4);
    x5 = param(5);
    C1 = exp(-x1);
    C2 = exp(-x2);
    C3 = exp(-x3);
    C4 = x4;
    C5 = x5;
    
    for i = 1: m
        res(i) = C5*C1^(C2^((1-C3^data(i))^C4))-y(i);
        J(i, 1) = -x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4) - 1)*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x1);
        J(i, 2) = -x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4 - 1)*exp(-x2)*log(exp(-x1))*(1 - exp(-x3)^data(i))^x4;
        J(i, 3) = x4*x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x3)^(data(i) - 1)*exp(-x3)*log(exp(-x1))*log(exp(-x2))*data(i)*(1 - exp(-x3)^data(i))^(x4 - 1);
        
        J(i, 4) = x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*log(1 - exp(-x3)^data(i))*log(exp(-x1))*log(exp(-x2))*(1 - exp(-x3)^data(i))^x4;
        J(i, 5) = exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4));
       
%         tmp(1,1) = x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4) - 1)*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x1) + x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4) - 2)*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-2*x1)*(exp(-x2)^((1 - exp(-x3)^data(i))^x4) - 1);
% 
%         tmp(1,2) = x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4) - 1)*exp(-x2)^((1 - exp(-x3)^data(i))^x4 - 1)*exp(-x1)*exp(-x2)*(1 - exp(-x3)^data(i))^x4 + x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4) - 1)*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x2)^((1 - exp(-x3)^data(i))^x4 - 1)*exp(-x1)*exp(-x2)*log(exp(-x1))*(1 - exp(-x3)^data(i))^x4;
% 
%         tmp(1,3) = - x4*x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4) - 1)*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x3)^(data(i) - 1)*exp(-x1)*exp(-x3)*log(exp(-x2))*data(i)*(1 - exp(-x3)^data(i))^(x4 - 1) - x4*x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4) - 1)*exp(-x2)^(2*(1 - exp(-x3)^data(i))^x4)*exp(-x3)^(data(i) - 1)*exp(-x1)*exp(-x3)*log(exp(-x1))*log(exp(-x2))*data(i)*(1 - exp(-x3)^data(i))^(x4 - 1);
% 
%         tmp(1,4) = - x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4) - 1)*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x1)*log(1 - exp(-x3)^data(i))*log(exp(-x2))*(1 - exp(-x3)^data(i))^x4 - x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4) - 1)*exp(-x2)^(2*(1 - exp(-x3)^data(i))^x4)*exp(-x1)*log(1 - exp(-x3)^data(i))*log(exp(-x1))*log(exp(-x2))*(1 - exp(-x3)^data(i))^x4;
% 
%         tmp(1,5) = -exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4) - 1)*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x1);
% 
%         tmp(2,1) = x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4 - 1)*exp(-x2)*(1 - exp(-x3)^data(i))^x4 + x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4) - 1)*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x2)^((1 - exp(-x3)^data(i))^x4 - 1)*exp(-x1)*exp(-x2)*log(exp(-x1))*(1 - exp(-x3)^data(i))^x4;
% 
%         tmp(2,2) = x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^(2*(1 - exp(-x3)^data(i))^x4 - 2)*exp(-2*x2)*log(exp(-x1))^2*(1 - exp(-x3)^data(i))^(2*x4) + x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4 - 1)*exp(-x2)*log(exp(-x1))*(1 - exp(-x3)^data(i))^x4 + x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4 - 2)*exp(-2*x2)*log(exp(-x1))*((1 - exp(-x3)^data(i))^x4 - 1)*(1 - exp(-x3)^data(i))^x4;
% 
%         tmp(2,3) = - x4*x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4 - 1)*exp(-x3)^(data(i) - 1)*exp(-x2)*exp(-x3)*log(exp(-x1))*data(i)*(1 - exp(-x3)^data(i))^(x4 - 1) - x4*x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4 - 1)*exp(-x3)^(data(i) - 1)*exp(-x2)*exp(-x3)*log(exp(-x1))*log(exp(-x2))*data(i)*(1 - exp(-x3)^data(i))^x4*(1 - exp(-x3)^data(i))^(x4 - 1) - x4*x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x2)^((1 - exp(-x3)^data(i))^x4 - 1)*exp(-x3)^(data(i) - 1)*exp(-x2)*exp(-x3)*log(exp(-x1))^2*log(exp(-x2))*data(i)*(1 - exp(-x3)^data(i))^x4*(1 - exp(-x3)^data(i))^(x4 - 1);
% 
%         tmp(2,4) = - x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4 - 1)*exp(-x2)*log(1 - exp(-x3)^data(i))*log(exp(-x1))*(1 - exp(-x3)^data(i))^x4 - x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4 - 1)*exp(-x2)*log(1 - exp(-x3)^data(i))*log(exp(-x1))*log(exp(-x2))*(1 - exp(-x3)^data(i))^(2*x4) - x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x2)^((1 - exp(-x3)^data(i))^x4 - 1)*exp(-x2)*log(1 - exp(-x3)^data(i))*log(exp(-x1))^2*log(exp(-x2))*(1 - exp(-x3)^data(i))^(2*x4);
% 
%         tmp(2,5) = -exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4 - 1)*exp(-x2)*log(exp(-x1))*(1 - exp(-x3)^data(i))^x4;
% 
%         tmp(3,1) = - x4*x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x3)^(data(i) - 1)*exp(-x3)*log(exp(-x2))*data(i)*(1 - exp(-x3)^data(i))^(x4 - 1) - x4*x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4) - 1)*exp(-x2)^(2*(1 - exp(-x3)^data(i))^x4)*exp(-x3)^(data(i) - 1)*exp(-x1)*exp(-x3)*log(exp(-x1))*log(exp(-x2))*data(i)*(1 - exp(-x3)^data(i))^(x4 - 1);
% 
%         tmp(3,2) = - x4*x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x3)^(data(i) - 1)*exp(-x3)*log(exp(-x1))*data(i)*(1 - exp(-x3)^data(i))^(x4 - 1) - x4*x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4 - 1)*exp(-x3)^(data(i) - 1)*exp(-x2)*exp(-x3)*log(exp(-x1))*log(exp(-x2))*data(i)*(1 - exp(-x3)^data(i))^x4*(1 - exp(-x3)^data(i))^(x4 - 1) - x4*x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x2)^((1 - exp(-x3)^data(i))^x4 - 1)*exp(-x3)^(data(i) - 1)*exp(-x2)*exp(-x3)*log(exp(-x1))^2*log(exp(-x2))*data(i)*(1 - exp(-x3)^data(i))^x4*(1 - exp(-x3)^data(i))^(x4 - 1);
% 
%         tmp(3,3) = x4^2*x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^(2*(1 - exp(-x3)^data(i))^x4)*exp(-x3)^(2*data(i) - 2)*exp(-2*x3)*log(exp(-x1))^2*log(exp(-x2))^2*data(i)^2*(1 - exp(-x3)^data(i))^(2*x4 - 2) + x4^2*x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x3)^(2*data(i) - 2)*exp(-2*x3)*log(exp(-x1))*log(exp(-x2))^2*data(i)^2*(1 - exp(-x3)^data(i))^(2*x4 - 2) - x4*x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x3)^(data(i) - 1)*exp(-x3)*log(exp(-x1))*log(exp(-x2))*data(i)*(1 - exp(-x3)^data(i))^(x4 - 1) + x4*x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x3)^(2*data(i) - 2)*exp(-2*x3)*log(exp(-x1))*log(exp(-x2))*data(i)^2*(1 - exp(-x3)^data(i))^(x4 - 2)*(x4 - 1) - x4*x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x3)^(data(i) - 2)*exp(-2*x3)*log(exp(-x1))*log(exp(-x2))*data(i)*(data(i) - 1)*(1 - exp(-x3)^data(i))^(x4 - 1);
% 
%         tmp(3,4) = x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x3)^(data(i) - 1)*exp(-x3)*log(exp(-x1))*log(exp(-x2))*data(i)*(1 - exp(-x3)^data(i))^(x4 - 1) + x4*x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x3)^(data(i) - 1)*exp(-x3)*log(1 - exp(-x3)^data(i))*log(exp(-x1))*log(exp(-x2))*data(i)*(1 - exp(-x3)^data(i))^(x4 - 1) + x4*x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^(2*(1 - exp(-x3)^data(i))^x4)*exp(-x3)^(data(i) - 1)*exp(-x3)*log(1 - exp(-x3)^data(i))*log(exp(-x1))^2*log(exp(-x2))^2*data(i)*(1 - exp(-x3)^data(i))^x4*(1 - exp(-x3)^data(i))^(x4 - 1) + x4*x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x3)^(data(i) - 1)*exp(-x3)*log(1 - exp(-x3)^data(i))*log(exp(-x1))*log(exp(-x2))^2*data(i)*(1 - exp(-x3)^data(i))^x4*(1 - exp(-x3)^data(i))^(x4 - 1);
% 
%         tmp(3,5) = x4*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x3)^(data(i) - 1)*exp(-x3)*log(exp(-x1))*log(exp(-x2))*data(i)*(1 - exp(-x3)^data(i))^(x4 - 1);
% 
%         tmp(4,1) = - x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*log(1 - exp(-x3)^data(i))*log(exp(-x2))*(1 - exp(-x3)^data(i))^x4 - x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4) - 1)*exp(-x2)^(2*(1 - exp(-x3)^data(i))^x4)*exp(-x1)*log(1 - exp(-x3)^data(i))*log(exp(-x1))*log(exp(-x2))*(1 - exp(-x3)^data(i))^x4;
% 
%         tmp(4,2) = - x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*log(1 - exp(-x3)^data(i))*log(exp(-x1))*(1 - exp(-x3)^data(i))^x4 - x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4 - 1)*exp(-x2)*log(1 - exp(-x3)^data(i))*log(exp(-x1))*log(exp(-x2))*(1 - exp(-x3)^data(i))^(2*x4) - x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x2)^((1 - exp(-x3)^data(i))^x4 - 1)*exp(-x2)*log(1 - exp(-x3)^data(i))*log(exp(-x1))^2*log(exp(-x2))*(1 - exp(-x3)^data(i))^(2*x4);
% 
%         tmp(4,3) = x4*x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x3)^(data(i) - 1)*exp(-x3)*log(1 - exp(-x3)^data(i))*log(exp(-x1))*log(exp(-x2))*data(i)*(1 - exp(-x3)^data(i))^(x4 - 1) - (x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x3)^(data(i) - 1)*exp(-x3)*log(exp(-x1))*log(exp(-x2))*data(i)*(1 - exp(-x3)^data(i))^x4)/(exp(-x3)^data(i) - 1) + x4*x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^(2*(1 - exp(-x3)^data(i))^x4)*exp(-x3)^(data(i) - 1)*exp(-x3)*log(1 - exp(-x3)^data(i))*log(exp(-x1))^2*log(exp(-x2))^2*data(i)*(1 - exp(-x3)^data(i))^x4*(1 - exp(-x3)^data(i))^(x4 - 1) + x4*x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x3)^(data(i) - 1)*exp(-x3)*log(1 - exp(-x3)^data(i))*log(exp(-x1))*log(exp(-x2))^2*data(i)*(1 - exp(-x3)^data(i))^x4*(1 - exp(-x3)^data(i))^(x4 - 1);
% 
%         tmp(4,4) = x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^(2*(1 - exp(-x3)^data(i))^x4)*log(1 - exp(-x3)^data(i))^2*log(exp(-x1))^2*log(exp(-x2))^2*(1 - exp(-x3)^data(i))^(2*x4) + x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*log(1 - exp(-x3)^data(i))^2*log(exp(-x1))*log(exp(-x2))^2*(1 - exp(-x3)^data(i))^(2*x4) + x5*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*log(1 - exp(-x3)^data(i))^2*log(exp(-x1))*log(exp(-x2))*(1 - exp(-x3)^data(i))^x4;
% 
%         tmp(4,5) = exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*log(1 - exp(-x3)^data(i))*log(exp(-x1))*log(exp(-x2))*(1 - exp(-x3)^data(i))^x4;
% 
%         tmp(5,1) = -exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4) - 1)*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x1);
% 
%         tmp(5,2) = -exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4 - 1)*exp(-x2)*log(exp(-x1))*(1 - exp(-x3)^data(i))^x4;
% 
%         tmp(5,3) = x4*exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*exp(-x3)^(data(i) - 1)*exp(-x3)*log(exp(-x1))*log(exp(-x2))*data(i)*(1 - exp(-x3)^data(i))^(x4 - 1);
% 
%         tmp(5,4) = exp(-x1)^(exp(-x2)^((1 - exp(-x3)^data(i))^x4))*exp(-x2)^((1 - exp(-x3)^data(i))^x4)*log(1 - exp(-x3)^data(i))*log(exp(-x1))*log(exp(-x2))*(1 - exp(-x3)^data(i))^x4;
% 
%         tmp(5,5) = 0;
%         H = H + res(i)*tmp;
    end
%     J = -J;
%     H = H + J'*J;
end