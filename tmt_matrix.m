function tm = tmt_matrix(i, beta, k0, ws, ns, pol)
    g2 = sqrt(beta ^ 2 - k0 ^ 2 * ns(i) ^ 2);
    g1 = sqrt(beta ^ 2 - k0 ^ 2 * ns(i - 1) ^ 2);
    theta = g2 * ws(i);
    rho = g2 / g1;
    if pol == "tm"
        rho = rho * ns(i - 1) ^ 2 / ns(i) ^ 2;
    end
    
    tm = .5 * [[(1 + rho) * exp(-theta), (1 - rho) * exp(theta)]; ...
               [(1 - rho) * exp(-theta), (1 + rho) * exp(theta)]] ;    
end