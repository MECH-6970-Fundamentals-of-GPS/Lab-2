%% PART I
clear; close all; clc;

addpath(genpath("./"));

data = load("RCVR_S1_data.mat").RCVR_S1;
ephem = data.ephem;
psrL1 = data.measurements.L1.psr;
time = data.GPS_time.seconds;

for i = 1:length(time)
    % Calculate current SV states
    [svPos, svVel, svB, svD, svPrns] = sv_positions(ephem, psrL1(i,:), time(i));
    sv_amt = 
    svPosLLA(:,:,i) = ecef2lla(svPos(:,:,i));
end

