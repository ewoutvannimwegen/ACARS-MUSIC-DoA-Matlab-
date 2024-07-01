%%
clear all;
close all;
% Uniform Linear Array of M isotropic antennas 
M = 4;
ula = phased.ULA('NumElements',M,'ElementSpacing',0.5);
%%
% Azimuth & Elevation angle of arrival for each Tx source antenna
ang1 = [40; 0];
ang2 = [-20; 0];
angs = [ang1 ang2];
c = physconst('LightSpeed');
fc = 300e6; 
lambda = c / fc;
pos = getElementPosition(ula)/lambda;
Nsamp = 1e3;
noisePwr = 0.01;
rs = rng(2007); 

% Simulate the received signal on the Rx antenna beam
sig = sensorsig(pos, Nsamp, angs, noisePwr);
%%
% Create your own MUSIC algorithm here!



%%
% Use Matlab's build-in MUSIC algo for comparison
music = phased.MUSICEstimator("SensorArray", ula, ...
    "OperatingFrequency", fc);

% Extract the angles of the Tx sources
ang = music(sig);
[~, ang1_est_idx] = max(ang);
ang(ang1_est_idx) = -Inf;
[~, ang2_est_idx] = max(ang);
est_angs = [ang1_est_idx, ang2_est_idx];
music.ScanAngles(est_angs)

% Plot the angular spectrum
plotSpectrum(music);
%%