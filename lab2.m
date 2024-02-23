%% PART I
clear; close all; clc;

% Constants
c = physconst('LightSpeed');
fL1 = 1575.42e6;
filename = 'RCVR_S1_data.mat';
toomersLLA = [32.606389, -85.481667, 214.65];
C = ECEF_ENU(toomersLLA(1), toomersLLA(2));
gamma = (77/60)^2;

% Load Data
addpath(genpath("./"));
data = load(filename).RCVR_S1;
ephem = data.ephem;
psrL1 = data.measurements.L1.psr;
dopL1 = data.measurements.L1.doppler;
psrL2 = data.measurements.L2.psr;
time = data.GPS_time.seconds;

X0 = zeros(8,1);
XuL1 = zeros(length(X0),length(time));
XuIF = zeros(length(X0),length(time));
[svPos, svVel, svB, svD, svPrns] = sv_positions(ephem, psrL1(1,:), time(1));
Xs = [svPos svB svVel svD];
rhoL1 = psrL1(1, svPrns)' + c*svB;
psr = (psrL2(1,:) - gamma*psrL1(1,:))./(1-gamma);
rhoIF = psr(svPrns)' + c*svB;
rho_dot = (-c/fL1)*dopL1(1, svPrns)' + c*svD;
XuL1(:,1) = gnssPVT(Xs, X0, rhoL1', rho_dot');
XuIF(:,1) = gnssPVT(Xs, X0, rhoIF', rho_dot');
for i = 2:length(time)
    [svPos, svVel, svB, svD, svPrns] = sv_positions(ephem, psrL1(i,:), time(i));
    Xs = [svPos svB svVel svD];
    rhoL1 = psrL1(i, svPrns)' + c*svB;
    psr = (psrL2(i,:) - gamma*psrL1(i,:))./(1-gamma);
    rhoIF = psr(svPrns)' + c*svB;
    rho_dot = (-c/fL1)*dopL1(i, svPrns)' + c*svD;
    [XuL1(:,i), HL1] = gnssPVT(Xs, XuL1(:,i-1), rhoL1', rho_dot');
    [XuIF(:,i), HIF] = gnssPVT(Xs, XuIF(:,i-1), rhoIF', rho_dot');
    uPosLLA = ecef2lla(XuL1(1:3,:)');
    DOP(i) = gpsStats(HL1, uPosLLA(1), uPosLLA(2));
end

uPosLLA = ecef2lla(XuL1(1:3,:)');
figure()
geoplot(uPosLLA(:,1),uPosLLA(:,2),'.')
geobasemap satellite

ENU = C*XuL1(1:3,:);
ENU_vel = C*XuL1(5:7,:);
pos_std = std(ENU,0,2);
vel_std = std(ENU_vel,0,2);

fprintf('Position Accuracy(in ENU): [%0.3g %0.3g %0.3g]m\n', pos_std);
fprintf('Velocity Accuracy(in ENU): [%0.3g %0.3g %0.3g]m\n', vel_std);

figure();
tiledlayout(3,1);
nexttile();
plot(time, ENU(1,:));
ylabel('East (m)');
ax = gca;
ax.FontSize = 16;
nexttile();
plot(time, ENU(2,:));
ylabel('North (m)');
ax = gca;
ax.FontSize = 16;
nexttile();
plot(time, ENU(3,:));
xlabel('Time (s)');
ylabel('Up (m)');
ax = gca;
ax.FontSize = 16;

figure();
tiledlayout(3,1);
nexttile();
plot(time, ENU_vel(1,:));
ylabel('East (m)');
ax = gca;
ax.FontSize = 16;
nexttile();
plot(time, ENU_vel(2,:));
ylabel('North (m)');
ax = gca;
ax.FontSize = 16;
nexttile();
plot(time, ENU_vel(3,:));
xlabel('Time (s)');
ylabel('Up (m)');
ax = gca;
ax.FontSize = 16;

%% PART II
clear; clc;

% Constants
c = physconst('LightSpeed');
fL1 = 1575.42e6;
filename = 'RCVR_D1_data.mat';
toomersLLA = [32.606389, -85.481667, 214.65];
C = ECEF_ENU(toomersLLA(1), toomersLLA(2));
gamma = (77/60)^2;

% Load Data
addpath(genpath("./"));
data = load(filename).RCVR_D1;
ephem = data.ephem;
psrL1 = data.measurements.L1.psr;
dopL1 = data.measurements.L1.doppler;
psrL2 = data.measurements.L2.psr;
dopL2 = data.measurements.L2.doppler;
time = data.GPS_time.seconds;

X0 = zeros(8,1);
XuL1 = zeros(length(X0),length(time));
XuIF = zeros(length(X0),length(time));
[svPos, svVel, svB, svD, svPrns] = sv_positions(ephem, psrL1(1,:), time(1));
Xs = [svPos svB svVel svD];
rhoL1 = psrL1(1, svPrns)' + c*svB;
psr = (psrL2(1,:) - gamma*psrL1(1,:))./(1-gamma);
rhoIF = psr(svPrns)' + c*svB;
rho_dot = (-c/fL1)*dopL1(1, svPrns)' + c*svD;
XuL1(:,1) = gnssPVT(Xs, X0, rhoL1', rho_dot');
XuIF(:,1) = gnssPVT(Xs, X0, rhoIF', rho_dot');
for i = 2:length(time)
    [svPos, svVel, svB, svD, svPrns] = sv_positions(ephem, psrL1(i,:), time(i));
    Xs = [svPos svB svVel svD];
    rhoL1 = psrL1(i, svPrns)' + c*svB;
    psr = (psrL2(i,:) - gamma*psrL1(i,:))./(1-gamma);
    rhoIF = psr(svPrns)' + c*svB;
    rho_dot = (-c/fL1)*dopL1(i, svPrns)' + c*svD;
    [XuL1(:,i), HL1] = gnssPVT(Xs, XuL1(:,i-1), rhoL1', rho_dot');
    [XuIF(:,i), HIF] = gnssPVT(Xs, XuIF(:,i-1), rhoIF', rho_dot');
    uPosLLA = ecef2lla(XuL1(1:3,:)');
    DOP(i) = gpsStats(HL1, uPosLLA(1), uPosLLA(2));
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