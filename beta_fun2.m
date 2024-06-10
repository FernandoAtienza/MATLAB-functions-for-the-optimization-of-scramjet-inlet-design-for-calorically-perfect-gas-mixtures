% Pressure Ratio Function

function f = beta_fun2(g, M, PR_guess)
    f = @(beta) - PR_guess + A(g, M, beta)^(g/(g-1)) * B(g, M, beta)^(1/(g-1));
end