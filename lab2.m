%% PART I
clear; close all; clc;

addpath(genpath("./"));

data = load("RCVR_S1_data.mat").RCVR_S1;
ephem = data.ephem;
psrL1 = data.measurements.L1.psr;
time = data.GPS_time.seconds;
c=physconst('LightSpeed');
fL1=1575.42e6;

for i = 1:length(time)
    % Part A: Calculate the sv pos as a function of GPS time
    [svPos, svVel, svB, svD, svPrns] = sv_positions(ephem, psrL1(i,:), time(i));

    % Part B: Calculate User Position
    psr=(data.measurements.L1.psr(i,svPrns)+c*svB')';
    doppler=data.measurements.L1.doppler(i,svPrns);
    psrR=((-c.*doppler./fL1)+c*svD')';
    [xhat(i,:),DOP(:,:,i)]=gpsPosVel(psr,[0,0,0,0],svPos,psrR,svVel, eye(length(svPrns)));
end

% Part C: Convert ECEF to LLA and Plot
uPosLLA=ecef2lla(xhat(:,1:3));
figure()
geoplot(uPosLLA(:,1),uPosLLA(:,2),'.')
geobasemap satellite

% Part D: Convert ECEF to ENU
time=time-time(1);
ToomersLLA=[32.606389,-85.481667,214.65];
wgs=wgs84Ellipsoid("meter");
[uPosE,uPosN,uPosU]=ecef2enu(xhat(:,1),xhat(:,2),xhat(:,3),ToomersLLA(1),ToomersLLA(2),ToomersLLA(3),wgs);
figure()
plot(time,uPosE)
title('East','.')
figure()
plot(time,uPosN)
title('North','.')
figure()
plot(time,uPosU)
title('Up','.')


% Part E: Calculate and find User Velocity in ENU
% Note: ECEF Velocity states computed in gpsPosVel function
[uVelE,uVelN,uVelU]=ecef2enu(xhat(:,4),xhat(:,5),xhat(:,6),ToomersLLA(1),ToomersLLA(2),ToomersLLA(3),wgs);





