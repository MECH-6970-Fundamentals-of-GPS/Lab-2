function [Xs, DOP] = gpsPosVel2(psr, X0, sPos, psrR, sVel, R)
% Performs Weighted Least Squares on GPS pseudoranges and range rates
%   Given satellite positions, velocites, and pseudoranges/range rates this 
%   function calculates the receiver position and velocity in the ECEF 
%   frame.
% Inputs:
%    psr  : pseudorange measurements [m,1]
%    uPos : initial user position in ECEF including bias [x,y,z,b] 
%               Note: if nothing is know about user location init to 
%               [0,0,0,0]
%    sPos : matrix of sv positions [m,3] - [x_sv,y_sv,z_sv]   
%    psrR : pseudorange rate measurements [m,1]
%               Note: This function assumes you have calculated psrR from
%               doppler and the center frequency of the band
%    R    : Measurement Noise Covariance Matrix [m,m]
%
% Outputs:
%    UserStates : The estimated user state [x,y,z,vx,vy,vz,b,d] 
%    DOP        : Dilution of Precision Matrix [8,8]

dX = 1e3 * ones(length(X0),1);
for i=1:100
    dr = sPos - X0(i,1:3);                % True Range
    r = vecnorm(dr, 2, 2);                  % Range Unit Vectors
    U = dr ./ r;                            % LOS Vector
    H = [-U, ones(length(psr),1)];          % Geometry Matrix
    rho_hat = r + X0(i,4);                % Estimated psr
    delpsr = psr - rho_hat;                 % Meas-Est Residual
    dX = (H'*R*H)\H'*R*delpsr;             % Least Squares
    X0(i+1,:)=X0(i,:)+dX';             % Newton-Raphson Update
    if (i == 100)
        disp('Solution did not converge\n')
        break
    end
    if norm(dX(1:3))<=1e-5
        break
    end
end
rrTilde=psrR-dot(sVel,U,2);                                     % Range Rate Meas
vel=(H'*H)\H'*rrTilde;                                          % Least Squares

Xs=zeros(1,8);
Xs(1,1:3) = X0(end,1:3);
Xs(1,4:6) = vel(1:3);
Xs(1,7) = X0(end,4);
Xs(1,8) = vel(4);
DOP=inv(H'*R*H);
end