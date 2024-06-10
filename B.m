% Function declarataion

% Function Coefficient

function value = B(g, M, beta)
    value = (g+1)/(2*g*((M * sind(beta))^2) - (g-1));
end