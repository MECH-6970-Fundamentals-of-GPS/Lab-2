%% PART I
clear; close all; clc;

% Constants
c = physconst('LightSpeed');
fL1 = 1575.42e6;
filename = 'RCVR_S1_data.mat';
toomersLLA = [32.606389, -85.481667, 214.65];
truthECEF = [422596.629, -5362864.287, 3415493.797];
truthLLA = ecef2lla(truthECEF);
Ctoomers = ECEF_ENU(toomersLLA(1), toomersLLA(2));
Ctruth = ECEF_ENU(truthLLA(1), truthLLA(2));
gamma = (77/60)^2;

% Load Data
addpath(genpath("./"));
data = load(filename).RCVR_S1;
ephem = data.ephem;
psrL1 = data.measurements.L1.psr;
dopL1 = data.measurements.L1.doppler;
psrVarL1 = data.measurements.L1.psr_variance;
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

    P_psr = diag(psrVarL1(i,svPrns))';
    P_pvt = inv(HL1'*inv(P_psr)*HL1);
    p_error_hat(i) = sqrt(sum(diag(P_pvt(1:3,1:3))));
end

uPosLLA = ecef2lla(XuL1(1:3,:)');
errorLLA = abs(uPosLLA - truthLLA);

figure()
geoplot(uPosLLA(:,1),uPosLLA(:,2),'*','LineWidth',3)
hold on;
title('Static Data Set')
geoplot(truthLLA(1), truthLLA(2),'rx','LineWidth',5,'MarkerSize',25);
geobasemap satellite
legend('L1 Only Estimate', 'Truth')
ax = gca;
ax.FontSize = 18;

ENU = Ctoomers*XuL1(1:3,:);
ENU_vel = Ctoomers*XuL1(5:7,:);
pos_std = std(ENU,0,2);
vel_std = std(ENU_vel,0,2);
time_std = std(XuL1(4,:));

fprintf('Position Accuracy(in ENU): [%0.3g %0.3g %0.3g]m\n', pos_std);
fprintf('Velocity Accuracy(in ENU): [%0.3g %0.3g %0.3g]m\n', vel_std);

figure();
plot(time, p_error_hat,'-x','MarkerSize',10);
xlabel('Time (s)');
ylabel('Position Error (m)');
title('Position Error Estimate vs. Time');
ax = gca;
ax.FontSize = 18;

figure();
plot(time, vecnorm(errorLLA'),'-x','MarkerSize',10);
xlabel('Time (s)');
ylabel('Position Error (m)');
title('Actual Position Error vs. Time');
ax = gca;
ax.FontSize = 18;

figure();
tiledlayout(3,1);
nexttile();
plot(time, ENU(1,:),'-x');
title('ENU Velocity v. Time');
subtitle('Origin at Toomers Corner, Auburn, AL');
ylabel('East (m)');
ax = gca;
ax.FontSize = 18;
nexttile();
plot(time, ENU(2,:),'-x');
ylabel('North (m)');
ax = gca;
ax.FontSize = 18;
nexttile();
plot(time, ENU(3,:),'-x');
xlabel('Time (s)');
ylabel('Up (m)');
ax = gca;
ax.FontSize = 18;

figure();
tiledlayout(3,1);
nexttile();
plot(time, ENU_vel(1,:),'-x');
title('ENU Velocity v. Time');
subtitle('Origin at Toomers Corner, Auburn, AL');
ylabel('East (m)');
ax = gca;
ax.FontSize = 18;
nexttile();
plot(time, ENU_vel(2,:),'-x');
ylabel('North (m)');
ax = gca;
ax.FontSize = 18;
nexttile();
plot(time, ENU_vel(3,:),'-x');
xlabel('Time (s)');
ylabel('Up (m)');
ax = gca;
ax.FontSize = 18;

%% PART II
clear;

% Constants
c = physconst('LightSpeed');
fL1 = 1575.42e6;
filename = 'RCVR_D1_data.mat';
toomersLLA = [32.606389, -85.481667, 214.65];
Ctoomers = ECEF_ENU(toomersLLA(1), toomersLLA(2));
gamma = (77/60)^2;

% Load Data
addpath(genpath("./"));
data = load(filename).RCVR_D1;

truthLLA = [data.true_pos.lat;data.true_pos.lon;data.true_pos.alt]';

ephem = data.ephem;
psrL1 = data.measurements.L1.psr;
dopL1 = data.measurements.L1.doppler;
psrVarL1 = data.measurements.L1.psr_variance;
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

    P_psr = diag(psrVarL1(i,svPrns))';
    P_pvt = inv(HL1'*inv(P_psr)*HL1);
    p_error_hat(i) = sqrt(sum(diag(P_pvt(1:3,1:3))));
end

uPosLLA = ecef2lla(XuL1(1:3,:)');
figure()
geoplot(truthLLA(:,1),truthLLA(:,2),'x','LineWidth',2,'MarkerSize',10);
hold('on');
geoplot(uPosLLA(:,1),uPosLLA(:,2),'*','LineWidth',2,'MarkerSize',10);
geobasemap satellite
legend('Truth','L1 Only')
ax = gca;
ax.FontSize = 18;

vel = Ctoomers*XuL1(5:7,:);
course = atan2(vel(1,:), vel(2,:));
avg_velocity = mean(vecnorm(vel));

fprintf('The average speed of the moving platform: %0.3g m/s\n', avg_velocity);

p_error = abs(truthLLA - uPosLLA);
figure();
plot(time, vecnorm(p_error'),'-x', 'MarkerSize',10);
xlabel('Time (s)');
ylabel('Position Error (m)');
title('True Position Error vs. Time');
subtitle('As Compared to Provided Truth');
ax = gca;
ax.FontSize = 18;

figure();
plot(time, p_error_hat,'-x','MarkerSize',10);
xlabel('Time (s)');
ylabel('Position Error (m)');
title('Position Error Estimate vs. Time');
ax = gca;
ax.FontSize = 18;

figure();
plot(time, vecnorm(vel),'-x');
title('Velocity vs. Time');
xlabel('Time (s)');
ylabel('Velocity m/s');
ax = gca;
ax.FontSize = 18;

figure();
plot(course,'-x');
title('GPS Course v. Time');
xlabel('Time (s)');
ylabel('Course(radians)');
ax = gca;
ax.FontSize = 18;

figure();
plot(XuL1(4,:),'-x');
title('Clock Bias vs. Time');
xlabel('Time (s)');
ylabel('Receiver Clock Bias (m)');
ax = gca;
ax.FontSize = 18;

figure();
plot(XuL1(8,:),'-x');
title('Clock Drift vs. Time');
xlabel('Time (s)');
ylabel('Receiver Clock Drift (m/s)');
ax = gca;
ax.FontSize = 18;

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
    DOP.H = H;
end

function [C] = ENU_ECEF(lat, lon)
    C = [-sind(lon), -cosd(lon)*sind(lat), cosd(lon)*cosd(lat);
          cosd(lon), -sind(lon)*sind(lat), sind(lon)*cosd(lat);
                  0,            cosd(lat),           sind(lat)];
end

function [C] = ECEF_ENU(lat, lon)
    C = ENU_ECEF(lat, lon)';
end