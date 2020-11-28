function quadratic(A)
    [x1, x2] = meshgrid(-5:0.2:5);
    z = 1/2*(x1.*(x1.*A(1,1)+x2.*A(2,1))+...
        x2.*(x1.*A(1,2)+x2.*A(2,2)));
    figure;
    subplot(122); 
    surf(x1, x2, z); 
    xlabel('x(1)'); ylabel('x(2)'); zlabel('1/2X^TAX');
    grid on;
    subplot(121); 
    contour(x1, x2, z, 10, 'LineWidth',1.2,'ShowText','on'); 
    xlabel('x(1)'); ylabel('x(2)');
    grid on;
end