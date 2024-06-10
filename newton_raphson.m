% Newton Raphson root finding method
%
% @author: Fernando Atienza González
%          Senior Undergraduate Aerospace Engineer Student
%          Universidad Carlos III de Madrid
%
% Last updated: 7th June 2024
% -------------------------------------------------------------------------
function [x,a] = newton_raphson(x0, f, df)

    a = 0;

    % Definitions (default values)
    it_max = 1000;
    tol = 1e-4;
    
    % Initialization
    it = 0;
    STOP = 1;

    % Loop
    while STOP > tol && it < it_max
        it = it + 1;
        
        % Evaluate f and df
        f0 = f(x0);
        df0 = df(x0);

        % Compute new value
        x = x0 - 2 * f(x0) / df(x0);
        
        % Re-estimation of first derivative
        f0_2 = f(x);

        % Compute solution
        x = x0 - f0^2 / (df0 * (f0 - f0_2));
        
        % Compute STOP
        STOP = max( abs((x - x0) / x), abs(f(x0) / x0) );

        % Update iteration
        x0 = x;
    end
    
    if STOP > tol || it > it_max
        a = 0;
        fprintf('***********************************************************\n')
        fprintf('Root algorithm not converged \n')
        fprintf('   Error       =  %8.2f [%%]  \n', STOP*100)
        fprintf('   Beta        =  %8.2f [deg]  \n', x)
        fprintf('   Iterations  =  %8.d [it] \n', it)
        fprintf('***********************************************************\n')
    end

end