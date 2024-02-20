function [UserStates,DOP] = gpsPosVel(psr,uPos,sPos,psrR,sVel,R)
%  Performs Weighted Least Squares on GPS pseudoranges and range rates
%   Given satellite positions, velocites, and pseudoranges this function
%   calculates the receiver position and velocity.
% Inputs:
%    psr  : pseudorange measurements [n,1]
%    uPos : initial user position including bias [x,y,z,b] 
%               Note: if nothing is know about user location init to 
%               [0,0,0,0]
%    sPos : matrix of sv positions [n,3] - [x_sv,y_sv,z_sv]   
%    psrR : pseudorange rate measurements [n,1]
%               Note: This function assumes you have calculated psrR from
%               doppler and the center frequency of the band
%    R    : Measurement Noise Covariance Matrix [n,n]
%
% Outputs:
%    UserStates : The estimated user state [x,y,z,vx,vy,vz,b,d]
%    DOP        : Dilution of Precision Matrix [n,n]

for i=1:100
    rdiff=sPos-uPos(i,1:3);                                     % True Range
    rhat=vecnorm(rdiff,2,2);                                    % Range Unit Vectors
    U=[rdiff(:,1)./rhat, rdiff(:,2)./rhat, rdiff(:,3)./rhat];   % LOS Vector
    H=[-U,ones(length(psr),1)];                                 % Geometry Matrix
    psrEst=rhat+uPos(i,4);                                      % Estimated psr
    delpsr=psr-psrEst;                                          % Meas-Est Error
    del=(H'*R*H)\H'*R*delpsr;                                   % Least Squares
    uPos(i+1,:)=uPos(i,:)+del';                                 % Newton- Raphson Update
    if (i == 100)
        disp('Solution did not converge\n')
        break
    end
    if norm(del(1:3))<=1e-5
        break
    end
end
rrTilde=psrR-dot(sVel,U,2);
vel=(H'*H)\H'*rrTilde;

UserStates=zeros(1,8);
UserStates(1,1:3) = uPos(end,1:3);
UserStates(1,4:6) = vel(1:3);
UserStates(1,7) = uPos(end,4);
UserStates(1,8) = vel(4);
DOP=inv(H'*R*H);




end