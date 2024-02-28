%% PART I
clear; close all; clc;

% Constants
c = physconst('LightSpeed');
fL1 = 1575.42e6;
lambdaL1 = 0.19;
filename = 'RCVR_S1_data.mat';
toomersLLA = [32.606389, -85.481667, 214.65];
staticPos = [422596.629, -5362864.287, 3415493.797];
staticPosLLA = ecef2lla(staticPos);
C = ECEF_ENU(toomersLLA(1), toomersLLA(2));
gamma = (77/60)^2;

% Load Data
addpath(genpath("./"));
data = load(filename).RCVR_S1;
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

% Setup of Vectors
X0 = zeros(8,1);
XuL1 = zeros(length(X0),length(time));
XuIF = zeros(length(X0),length(time));
XuW = zeros(length(X0),length(time));
XuPD = zeros(length(X0),length(time));
XuBEST = zeros(length(X0),length(time));
[svPos, svVel, svB, svD, svPrns] = sv_positions(ephem, psrL1(1,:), time(1));
Xs = [svPos svB svVel svD];

% Normal L1 Measurements
rhoL1 = psrL1(1, svPrns)' + c*svB;
dop = (-c/fL1)*dopL1(1, :)';
rho_dot = dop(svPrns) + c*svD;

% Dual Frequency Ionospheric Correction
psr = (psrL2(1,:) - gamma*psrL1(1,:))./(1-gamma);
rhoL12 = psrL1(2, svPrns)' + c*svB;
rhoIF = psr(svPrns)' + c*svB;

% Pseudorange Doppler Calculation
psrdop = (psrL1(2,:) - psrL1(1,:))';

rho_dotc = dop;
rho_dotc(~isnan(psrdop)) = psrdop(~isnan(psrdop));
rho_dotc = rho_dotc(svPrns) + c*svD;

% Initial Position Estimate
XuL1(:,1) = gnssPVT(Xs, X0, rhoL1', rho_dot');
XuIF(:,1) = gnssPVT(Xs, X0, rhoIF', rho_dot');
XuW(:,1) = gnssPVT(Xs, X0,  rhoL1', rho_dot');
XuPD(:,1) = gnssPVT(Xs, X0, rhoL1', rho_dotc');
XuBEST(:,1) = gnssPVT(Xs, X0, rhoL1', rho_dotc');

for i = 2:length(time)
    [svPos, svVel, svB, svD, svPrns] = sv_positions(ephem, psrL1(i,:), time(i));
    Xs = [svPos svB svVel svD];

    % Normal L1 Measurements
    rhoL1 = psrL1(i, svPrns)' + c*svB;
    dop = (-c/fL1)*dopL1(i, :)';
    rho_dot = dop(svPrns) + c*svD;

    % Dual Frequency Ionospheric Correction
    psr(i,:) = (psrL2(i,:) - gamma*psrL1(i,:))./(1-gamma);
    rhoIF = psr(i, svPrns)' + c*svB;
    IFPrns = ~isnan(rhoIF);

    % Weighted Calculation
    r = (db2mag(cnoL1(i,svPrns))./psrVarL1(i,svPrns));
    R = diag(r);

    % Pseudorange Doppler Calculation
    psrdop = (psrL1(i,:) - psrL1(i-1,:))';

    rho_dotc = dop;
    rho_dotc(~isnan(psrdop)) = psrdop(~isnan(psrdop));
    rho_dotc = rho_dotc(svPrns) + c*svD;

    [XuL1(:,i), HL1] = gnssPVT(Xs, XuL1(:,i-1), rhoL1', rho_dot');
    [XuIF(:,i), HIF] = gnssPVT(Xs(IFPrns,:), XuIF(:,i-1), ...
        rhoIF(IFPrns)', rho_dot(IFPrns)');
    [XuW(:,i), HW] = gnssPVT(Xs, XuL1(:,i-1), rhoL1', rho_dot',R);
    [XuPD(:,i), HCD] = gnssPVT(Xs, XuL1(:,i-1), rhoL1', rho_dotc');
    [XuBEST(:,i), HBEST] = gnssPVT(Xs, XuL1(:,i-1), rhoL1', rho_dotc',R);
end

uPosLLA = ecef2lla(XuL1(1:3,:)');
uPosLLAIF = ecef2lla(XuIF(1:3,:)');
uPosLLAW = ecef2lla(XuW(1:3,:)');
figure()
geoplot(uPosLLA(:,1),uPosLLA(:,2),'x','LineWidth',3);
hold on
title('Static Data Set')
geoplot(uPosLLAIF(:,1),uPosLLAIF(:,2),'o','LineWidth',3);
geoplot(uPosLLAW(:,1),uPosLLAW(:,2),'*','LineWidth',3);
geoplot(staticPosLLA(1),staticPosLLA(2),'xb','LineWidth',10,'MarkerSize',50)
geobasemap satellite
legend('L1 Only', 'L1/L2', 'Weighted', 'Truth');
ax = gca;
ax.FontSize = 18

pos_std = std(C*XuL1(1:3,:),0,2);
vel_std = std(C*XuL1(5:7,:),0,2);
pos_stdIF = std(C*XuIF(1:3,:),0,2);
vel_stdIF = std(C*XuIF(5:7,:),0,2);
pos_stdW = std(C*XuW(1:3,:),0,2);
vel_stdW = std(C*XuW(5:7,:),0,2);
pos_stdPD = std(C*XuPD(1:3,:),0,2);
vel_stdPD = std(C*XuPD(5:7,:),0,2);
pos_stdBEST = std(C*XuBEST(1:3,:),0,2);
vel_stdBEST = std(C*XuBEST(5:7,:),0,2);

fprintf('L1 Only Position Accuracy(in ENU): [%0.3g %0.3g %0.3g]m\n', pos_std);
fprintf('L1 Only Velocity Accuracy(in ENU): [%0.3g %0.3g %0.3g]m\n', vel_std);
fprintf('L1/L2 Position Accuracy(in ENU): [%0.3g %0.3g %0.3g]m\n', pos_stdIF);
fprintf('L1/L2 Velocity Accuracy(in ENU): [%0.3g %0.3g %0.3g]m\n', vel_stdIF);
fprintf('Weighted L1 Only Position Accuracy(in ENU): [%0.3g %0.3g %0.3g]m\n', pos_stdW);
fprintf('Weighted L1 Only Velocity Accuracy(in ENU): [%0.3g %0.3g %0.3g]m\n', vel_stdW);
fprintf('L1 + Pseudorange Doppler Position Accuracy(in ENU): [%0.3g %0.3g %0.3g]m\n', pos_stdPD);
fprintf('L1 + Pseudorange Doppler Velocity Accuracy(in ENU): [%0.3g %0.3g %0.3g]m\n', vel_stdPD);
fprintf('Weighted L1 + Pseudorange Doppler Position Accuracy(in ENU): [%0.3g %0.3g %0.3g]m\n', pos_stdBEST);
fprintf('Weighted L1 + Pseudorange Doppler Velocity Accuracy(in ENU): [%0.3g %0.3g %0.3g]m\n', vel_stdBEST);

%% PART II
clear;

% Constants
c = physconst('LightSpeed');
fL1 = 1575.42e6;
lambdaL1 = 0.19;
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
time = data.GPS_time.seconds;

% Setup of Vectors
X0 = zeros(8,1);
XuL1 = zeros(length(X0),length(time));
XuW = zeros(length(X0),length(time));
XuPD = zeros(length(X0),length(time));
XuBEST = zeros(length(X0),length(time));
[svPos, svVel, svB, svD, svPrns] = sv_positions(ephem, psrL1(1,:), time(1));
Xs = [svPos svB svVel svD];

% Normal L1 Measurements
rhoL1 = psrL1(1, svPrns)' + c*svB;
dop = (-c/fL1)*dopL1(1, :)';
rho_dot = dop(svPrns) + c*svD;

% Pseudorange Doppler Calculation
psrdop = (psrL1(2,:) - psrL1(1,:))';

rho_dotc = dop;
rho_dotc(~isnan(psrdop)) = psrdop(~isnan(psrdop));
rho_dotc = rho_dotc(svPrns) + c*svD;

% Initial Position Estimate
XuL1(:,1) = gnssPVT(Xs, X0, rhoL1', rho_dot');
XuW(:,1) = gnssPVT(Xs, X0,  rhoL1', rho_dot');
XuPD(:,1) = gnssPVT(Xs, X0, rhoL1', rho_dot');
XuBEST(:,1) = gnssPVT(Xs, X0, rhoL1', rho_dotc');

for i = 2:length(time)
    [svPos, svVel, svB, svD, svPrns] = sv_positions(ephem, psrL1(i,:), time(i));
    Xs = [svPos svB svVel svD];

    % Normal L1 Measurements
    rhoL1 = psrL1(i, svPrns)' + c*svB;
    dop = (-c/fL1)*dopL1(i, :)';
    rho_dot = dop(svPrns) + c*svD;

    % Weighted Calculation
    r = db2mag(cnoL1(i,svPrns));
    R = diag(r);

    % Pseudorange Doppler Calculation
    psrdop = (psrL1(i,:) - psrL1(i-1,:))';

    rho_dotc = dop;
    rho_dotc(~isnan(psrdop)) = psrdop(~isnan(psrdop));
    rho_dotc = rho_dotc(svPrns) + c*svD;

    [XuL1(:,i), HL1] = gnssPVT(Xs, XuL1(:,i-1), rhoL1', rho_dot');
    [XuW(:,i), HW] = gnssPVT(Xs, XuL1(:,i-1), rhoL1', rho_dot',R);
    [XuPD(:,i), HCD] = gnssPVT(Xs, XuL1(:,i-1), rhoL1', rho_dotc');
    [XuBEST(:,i), HBEST] = gnssPVT(Xs, XuL1(:,i-1), rhoL1', rho_dotc',R);
end

uPosLLA = ecef2lla(XuL1(1:3,:)');
uPosLLAW = ecef2lla(XuW(1:3,:)');
figure()
geoplot(uPosLLA(:,1),uPosLLA(:,2),'x','LineWidth',3, 'MarkerSize',25);
hold on
title('Dynamic Data Set')
geoplot(uPosLLAW(:,1),uPosLLAW(:,2),'*','LineWidth',3, 'MarkerSize',25);
geobasemap satellite
legend('L1 Only', 'Weighted');
ax = gca;
ax.FontSize = 18;

ENU = C*XuL1(1:3,:);
ENUvelL1 = C*XuL1(5:7,:);
ENUW = C*XuW(1:3,:);
ENUvelW = C*XuW(5:7,:);
ENUPD = C*XuPD(1:3,:);
ENUvelPD = C*XuPD(5:7,:);
ENUBEST = C*XuBEST(1:3,:);
ENUvelBEST = C*XuBEST(5:7,:);

avgVelL1 = mean(vecnorm(ENUvelL1));
avgVelW = mean(vecnorm(ENUvelW));
avgVelCD = mean(vecnorm(ENUvelPD));
avgVelBEST = mean(vecnorm(ENUvelBEST));

fprintf('L1 Only Average Velocity: %0.3g m/s\n', avgVelL1);
fprintf('Weighted L1 Average Velocity: %0.3g m/s\n', avgVelW);
fprintf('L1 + Psuedorange Doppler Average Velocity: %0.3g m/s\n', avgVelCD);
fprintf('Weighted L1 + Psuedorange Doppler Average Velocity: %0.3g m/s\n', avgVelBEST);

figure();
hold('on');
plot(time, vecnorm(ENUvelL1), '-x');
plot(time, vecnorm(ENUvelW), '-o');
plot(time, vecnorm(ENUvelPD), '-*');
plot(time, vecnorm(ENUvelBEST), '-^');
title('Velocity vs. Time');
xlabel('Time (s)');
ylabel('Velocity m/s');
ax = gca;
ax.FontSize = 16;
legend('L1 Only', 'Weighted', 'Pseudorange Doppler', ...
    'Weighted + Pseudorange Doppler');

figure();
hold('on');
plot(XuL1(4,:), '-x');
plot(XuW(4,:), '-o');
plot(XuPD(4,:), '-*');
plot(XuBEST(4,:), '-^');
title('Clock Bias vs. Time');
xlabel('Time (s)');
ylabel('Receiver Clock Bias (m)');
legend('L1 Only', 'Weighted', 'Pseudorange Doppler', ...
    'Weighted + Pseudorange Doppler');
ax = gca;
ax.FontSize = 16;

figure();
hold('on');
plot(XuL1(8,:), '-x');
plot(XuW(8,:), '-o');
plot(XuPD(8,:), '-*');
plot(XuBEST(8,:), '-^');
title('Clock Drift vs. Time');
xlabel('Time (s)');
ylabel('Receiver Clock Drift (m/s)');
legend('L1 Only', 'Weighted', 'Pseudorange Doppler', ...
    'Weighted + Pseudorange Doppler');
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