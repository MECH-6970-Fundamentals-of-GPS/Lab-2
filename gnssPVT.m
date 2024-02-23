function [Xu, H] = gnssPVT(Xs, X0, rho, rho_dot)
    Xu = X0;
    dX = 1e3 * ones(4,1);
    while norm(dX) > 1e-4
        dr = Xs(:,1:3) - Xu(1:3)';
        r = vecnorm(dr, 2, 2);
        rho_hat = r' + Xu(4);
        U = dr./r;
        H = [-U ones(length(U),1)];
        dX = pinv(H)*(rho - rho_hat)';
        Xu(1:4) = Xu(1:4) + dX;
    end
    r_dot = rho_dot - dot(Xs(:,5:7), U, 2)';
    Xu(5:8) = pinv(H)*r_dot';
end