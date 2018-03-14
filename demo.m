function demo_fminlev()
    close all;

    A = 0.1; delta = 0; w_d = 1; phi_0 = 0;             % not so good, still works
    % A = 2.47; delta = 0.12; w_d = 2.3; phi_0 = 5.13;  % pretty good initialization
    % A = 0; delta = 0; w_d = 0; phi_0 = 0;             % zero intiialization: fails !

    x0 = [A; delta; w_d; phi_0];

    X = zeros(1, 20); Y = zeros(1, 20);
    X(1, 1:10) = [ 0.23 0.81  1.28  1.88  2.35 2.96 3.53 3.68  4.12  4.80 ];
    Y(1, 1:10) = [ 2.33 2.53 -0.36 -2.80 -1.49 1.93 2.12 1.51 -0.78 -2.09 ];
    X(1, 11:20) = [  5.27 5.73 6.16 6.43 6.69  7.12  7.76 8.21 8.78 9.38 ];
    Y(1, 11:20) = [ -0.32 1.62 1.86 1.11 0.04 -1.47 -1.07 0.51 1.44 0.15 ];

    figure; scatter(X, Y); grid on; title('Levenberg Marquardt Optimization'); hold on;
    XX = 0:0.01:10;
    h = [];

    function E = fenergy(x)
        f = x(1)*exp(- x(2)*X).*cos(x(3)*X + x(4));
        E = 0.5*sum( (f - Y).^2 );
    end

    function [H,f,c] = flinear(x0)
        f0  = x0(1)*exp(- x0(2)*X).*cos(x0(3)*X + x0(4));
        df1 = exp(- x0(2)*X).*cos(x0(3)*X + x0(4));
        df2 = - x0(1)*exp(- x0(2)*X).*cos(x0(3)*X + x0(4)).*X;
        df3 = - x0(1)*exp(- x0(2)*X).*sin(x0(3)*X + x0(4)).*X;
        df4 = - x0(1)*exp(- x0(2)*X).*sin(x0(3)*X + x0(4));
        A = [df1; df2; df3; df4]'; b = (f0 - Y)';
        H = A'*A; f = A'*b; c = b'*b;
    end

    function fshowiter(x,E,mu,progress,iter)
        YY = x(1)*exp(- x(2)*XX).*cos(x(3)*XX + x(4));
        delete(h);
        h = plot(XX,YY,'color','r');
        axis([0 10 -3 4]); drawnow;
        fprintf('iter=%02i a=%s E=%f mu=%f progress= %f\n', iter, x2str(x,'%+4.2f'), E, mu, progress);
        pause(.3);
    end

    tolerance = 1e-5; iterations = 1000;
    fminlev(@fenergy, @flinear, @fshowiter, tolerance, iterations, x0);

end


function str = x2str(x,formatSpec)
    if nargin < 2, formatSpec = []; end
    if ischar(x)
        str = x;
    else if isvector(x) && length(x) > 1
            str = sprintf('(%s',num2str(x(1),formatSpec));
            for j = 2:length(x)
                str = strcat(str,[' ' num2str(x(j),formatSpec)]);
            end
            str = strcat(str,')');
        else
            str = num2str(x,formatSpec);
        end
    end
end
