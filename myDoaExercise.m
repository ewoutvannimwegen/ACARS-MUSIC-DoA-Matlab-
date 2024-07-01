%%
clear all
clc
close all

% Uniform Linear Array of M isotropic antennas 
c = 3e8;
fc = 100e6; 
lambda = c / fc
M = 4;
dsep=lambda;
ula = phased.ULA('NumElements',M,'ElementSpacing',dsep);
%%
% Azimuth & Elevation angle of arrival for each Tx source antenna
ang1 = [40; 0];
ang2 = [-20; 0];
ang3 = [-10; 0];
angs = [ang1 ang2];

pos = getElementPosition(ula)/lambda
Nsamp = 1e3;
noisePwr = 0.01;
rs = rng(2007); 

% Simulate the received signal on the Rx antenna beam
sig = sensorsig(pos, Nsamp, angs, noisePwr);
%%
% My implementation of MUSIC

mySig = sig';

% Covariance matrix Rx of all Rx sensors
% if M is the number of sensors
% the covariance matrix should be MxM
sigH = conj(mySig');
Rx = (mySig * sigH) / Nsamp;

% Eigenvalue decompostion
[eigVecs, eigValsDM] = eig(Rx);

% Extract only items on diagonal, sort and return index
[~, idx] = sort(diag(eigValsDM), 'descend');

% Sort the eigenvectors according to the eigenvalues
eigVecs = eigVecs(:, idx);

% Extract the noise subspace consisting of
% the eigenvectors corresponding to the smallest M - K eigenvalues
Un = eigVecs(:, end-M+1:end);
UnH = conj(Un');

% Search through all angles in search space
thetaScan = -90:0.1:90;
Pmusic = zeros(size(thetaScan));
for i = 1:length(thetaScan)
    k = (2*pi)/lambda;
    steerVec = exp(-1j * k * dsep * (0:M-1)' * sind(thetaScan(i)));
    steerVecH = conj(steerVec');
    Pmusic(i) = 1 / (steerVecH * Un * UnH * steerVec);
end

figure;
plot(thetaScan, 10*log(abs(Pmusic).^2));
xlabel("Angle of Arrival [deg]");
ylabel("Power [dB]");
title("MUSIC DoA Estimation");
grid on;
xticks(-90:10:90);
xlim([-90, 90]);
ylim([0, 100]);

% Find the peaks
[peaks, peaksIdxs] = findpeaks(abs(Pmusic));
[~, sortedIndices] = sort(peaks, 'descend');
sortedPeaks = peaks(sortedIndices)
sortedPeaksIdxs = peaksIdxs(sortedIndices);
DoA = thetaScan(sortedPeaksIdxs(1:length(angs)))

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
doaRef = music.ScanAngles(est_angs)

% Plot the angular spectrum
%plotSpectrum(music);
%%