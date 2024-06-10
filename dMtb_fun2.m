% 1st Derivative of M-θ-β relationship function

function df = dMtb_fun2(g, M)
    df = @(beta) 2 * (-4 * M^2 + M^4 * (1 + g * cosd(2 * beta) ) + (2 + M^2 * (1 + g) ) * cscd(beta)^2 ) / ...
        (2 + M^2 * g + M^2 * cosd(2 * beta))^2;
end