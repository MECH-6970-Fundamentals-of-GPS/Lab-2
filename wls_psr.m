function Xo = wls_psr(svPos, svVel, Zmeas, R)
%wls_psr.m Performs Weighted Least Squares on GPS pseudoranges
%   Given satellite positions, velocites, and pseudoranges this function
%   calculates the receiver position and velocity.
% Inputs:
%    svPos  : a matrix of satellite positions [x,y,z]
%    svVel  : a matrix of satellite velocities [vx,vy,vz]
%    Zmeas  : a vector of pseudoranges and doppler measurements
%    R      : Measurement Noise Covariance Matrix
%
% Outputs:
%    Xo     : The estimated user state [x,y,z,vx,vy,vz,b,d]

    Xo = zeros(8,1);                                    % User State
    % Iterate until X converges
    dX = ones(length(Xo),1)*100;                        % Change in User State
    while(norm(dX) > 1e-6)
        dr = svPos - Xo(1:3)';                          % True range (vector) [m]
        r = vecnorm(dr,2,2);                            % True range (scalar) [m]
        U = dr ./ r;                                    % Line of Sight Unit Vectors
        Zhat = [r + Xo(7);
                sum(U.*(svVel - Xo(4:6)'),2) + Xo(8)];  % Estimate pseudoranges [m]

        L = length(U);
        H = [-U, zeros(L,3), ones(L,1), zeros(L,1); ...
             zeros(L,3), -U, zeros(L,1), ones(L,1)];    % Geometry Matrix
        dZ = Zmeas - Zhat;                              % Difference in Measurement & Estimate
        dX = (H'*R*H)^(-1) * (H'*R*dZ);                 % State Update
        Xo = Xo + dX;                                   % Update State
    end
end