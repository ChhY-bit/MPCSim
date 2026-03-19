function terminal_term = MPCSim_tmnl(prob, Kt, alpha)
    At = prob.Ad+prob.Bd*Kt;
    n = size(At,1);
    St = dlyap(At',eye(n));
    [V,D]=eig(St); V = V';
    sigma = sum(D)';
    beta = sqrt(alpha/n)*(1./sqrt(sigma));
    P = dlyap(At',Kt'*prob.R_each*Kt+At'*prob.Q_each*At);

    terminal_term.P = P;
    terminal_term.beta = beta;
    terminal_term.V = V;
end