% M-θ-β relationship function

function f = Mtb_fun2(g, M, theta)
    f = @(beta) -tand(theta) + (2 * cotd(beta) * ((M * sind(beta))^2 - 1) ./ (((g + cosd(2 * beta)) * M^2) + 2));
end