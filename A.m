% Function Declarataion

% Function coefficient

function value = A(g, M, beta)
    value = ((g+1)*((M * sind(beta))^2))/((g-1)*((M * sind(beta))^2) +2);
end