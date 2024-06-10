% 1st Derivative of Pressure Ratio Function

function df = dbeta_fun2(g, M)
    df = @(beta)- (4 * g * A(g, M, beta)^( (2*g - 1) / (g - 1) ) * B(g, M, beta)^( (2 - g) / (g - 1) ) * ...
        cotd(beta) * ( cscd(beta) - M^2 * sind(beta) )^2 ) / ((2 * M^3 * g * sind(beta)^2 - M * (g - 1) )^2);
end