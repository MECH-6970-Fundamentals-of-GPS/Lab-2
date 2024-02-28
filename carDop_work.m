%% PART I
clear; close all; clc;

% Constants
c = physconst('LightSpeed');
fL1 = 1575.42e6;
lambdaL1 = c/fL1;
filename = 'RCVR_D1_data.mat';
toomersLLA = [32.606389, -85.481667, 214.65];
C = ECEF_ENU(toomersLLA(1), toomersLLA(2));
gamma = (77/60)^2;

% Load Data
addpath(genpath("./"));
data = load(filename).RCVR_D1;
ephem = data.ephem;
psrL1 = data.measurements.L1.psr;
psrVarL1 = data.measurements.L1.psr_variance;
dopL1 = data.measurements.L1.doppler;
carL1 = data.measurements.L1.carrier_phase;
cnoL1 = data.measurements.L1.carrier_to_noise;
psrL2 = data.measurements.L2.psr;
psrVarL2 = data.measurements.L2.psr_variance;
carL2 = data.measurements.L2.carrier_phase;
cnoL2 = data.measurements.L2.carrier_to_noise;
time = data.GPS_time.seconds;

X0 = zeros(8,1);
XuL1 = zeros(length(X0),length(time));
XuIF = zeros(length(X0),length(time));
[svPos, svVel, svB, svD, svPrns] = sv_positions(ephem, psrL1(1,:), time(1));
Xs = [svPos svB svVel svD];

rhoL1 = psrL1(1, svPrns)' + c*svB;
range_rate(1,:) = rhoL1;
rhoL12 = psrL1(2, svPrns)' + c*svB;
rho_dot = (-c/fL1)*dopL1(1, :)';

phiSD = -(carL1(2,:) - carL1(1,:))'./(2*pi);
psrSD = (psrL1(2,:) - psrL1(1,:))';
N = psrSD - phiSD;
rho_dotc = (phiSD + N);

rho_dot(~isnan(rho_dotc)) = rho_dotc(~isnan(rho_dotc));
rho_dot = rho_dot(svPrns) + c*svD;

XuL1(:,1) = gnssPVT(Xs, X0, rhoL1', rho_dot');
for i = 2:length(time)
    [svPos, svVel, svB, svD, svPrns] = sv_positions(ephem, psrL1(i,:), time(i));
    Xs = [svPos svB svVel svD];

    psr(i,:) = (psrL2(i,:) - gamma*psrL1(i,:))./(1-gamma);
    rhoL1 = psrL1(i, svPrns)' + c*svB;
    rho_dot = (-c/fL1)*dopL1(i, :)';

    phiSD = -(carL1(i,:) - carL1(i-1,:))'./(2*pi);
    psrSD = (psrL1(i,:) - psrL1(i-1,:))';
    N = psrSD - phiSD;
    rho_dotc = (phiSD + N);
    
    rho_dot(~isnan(rho_dotc)) = rho_dotc(~isnan(rho_dotc));
    rho_dot = rho_dot(svPrns) + c*svD;

    r = db2mag(cnoL1(i,svPrns));
    R = diag(r);

    % [XuL1(:,i), HL1] = gnssPVT(Xs(prns,:), XuL1(:,i-1), rhoL1(prns)', rho_dot');
    [XuL1(:,i), HL1] = gnssPVT(Xs, XuL1(:,i-1), rhoL1', rho_dot',R);
end

uPosLLA = ecef2lla(XuL1(1:3,:)');
figure()
geoplot(uPosLLA(:,1),uPosLLA(:,2),'.')
geobasemap satellite

vel = C*XuL1(5:7,:);
course = atan2(vel(1,:), vel(2,:));
avg_velocity = mean(vecnorm(vel));

fprintf('The average speed of the moving platform: %0.3g m/s\n', avg_velocity);

figure();
plot(time, vecnorm(vel));
title('Velocity vs. Time');
xlabel('Time (s)');
ylabel('Velocity m/s');
ax = gca;
ax.FontSize = 16;

figure();
plot(course);
title('GPS Course v. Time');
xlabel('Time (s)');
ylabel('Course(radians)');
ax = gca;
ax.FontSize = 16;

figure();
plot(XuL1(4,:));
title('Clock Bias vs. Time');
xlabel('Time (s)');
ylabel('Receiver Clock Bias (m)');
ax = gca;
ax.FontSize = 16;

figure();
plot(XuL1(8,:));
title('Clock Drift vs. Time');
xlabel('Time (s)');
ylabel('Receiver Clock Drift (m/s)');
ax = gca;
ax.FontSize = 16;


%% FUNCTIONS
function [DOP] = gpsStats(H, lat, lon)
    C = ECEF_ENU(lat,lon);
    HH = inv(H'*H);
    HH(1:3, 1:3) = C*HH(1:3, 1:3)*C';
    DOP.G = sqrt(sum(diag(HH)));
    DOP.P = sqrt(sum(diag(HH(1:3,1:3))));
    DOP.H = sqrt(sum(diag(HH(1:2,1:2))));
    DOP.V = sqrt(HH(3,3));
    DOP.T = sqrt(HH(4,4));
end

function [C] = ENU_ECEF(lat, lon)
    C = [-sind(lon), -cosd(lon)*sind(lat), cosd(lon)*cosd(lat);
          cosd(lon), -sind(lon)*sind(lat), sind(lon)*cosd(lat);
                  0,            cosd(lat),           sind(lat)];
end

function [C] = ECEF_ENU(lat, lon)
    C = ENU_ECEF(lat, lon)';
end