%% PART I
clear; close all; clc;

% Constants
c = physconst('LightSpeed');
fL1 = 1575.42e6;
filename = 'RCVR_S1_data.mat';
toomersLLA = [32.606389, -85.481667, 214.65];

% Load Data
addpath(genpath("./"));
data = load(filename).RCVR_S1;
ephem = data.ephem;
psrL1 = data.measurements.L1.psr;
dopL1 = data.measurements.L1.doppler;
psrL2 = data.measurements.L2.psr;
time = data.GPS_time.seconds;

X0 = [0 0 0 0];

for i = 1:length(time)
    % Part A: Calculate the sv pos as a function of GPS time
    [svPos, svVel, svB, svD, svPrns] = sv_positions(ephem, psrL1(i,:), time(i));

    % Part B: Calculate User Position
    rho = psrL1(i, svPrns)' + c*svB;
    rho_dot = -c.*dopL1(i, svPrns)'./fL1 + c*svD;
    [Xu(i,:), DOP(:,:,i)] = gpsPosVel2(rho, X0, svPos, rho_dot, svVel, eye(length(svPrns)));
end

% Part C: Convert ECEF to LLA and Plot
uPosLLA = ecef2lla(Xu(:,1:3));
figure()
geoplot(uPosLLA(:,1),uPosLLA(:,2),'.')
geobasemap satellite

% % Part D: Convert ECEF to ENU
% time = time - time(1);
% wgs = wgs84Ellipsoid("meter");
% [uPosE,uPosN,uPosU] = ecef2enu(Xu(:,1), Xu(:,2), Xu(:,3), toomersLLA(1), toomersLLA(2), toomersLLA(3),wgs);
% figure()
% plot(time,uPosE)
% title('East','.')
% figure()
% plot(time,uPosN)
% title('North','.')
% figure()
% plot(time,uPosU)
% title('Up','.')
% 
% % Part E: Calculate and find User Velocity in ENU
% % Note: ECEF Velocity states computed in gpsPosVel function
% [uVelE,uVelN,uVelU] = ecef2enu(Xu(:,4),Xu(:,5),Xu(:,6),toomersLLA(1),toomersLLA(2),toomersLLA(3),wgs);





