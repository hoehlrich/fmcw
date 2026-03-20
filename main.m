params.fs = 20e6;       % Sampling frequency
params.Tc = 40e-6;      % Chirp period
params.B = 0.4e9;       % Chirp bandwidth
params.fc = 77e9;       % Carrier frequency
params.Nc = 128;        % Number of chirps
params.noise = 0;       % Noise amplitude
params.swerling = true; % Swerling enabled


ranges = [50,60,120];
velocities = [5, -10, -15];
[ranges, velocities] = spread_targets(ranges, velocities, 0.25, 0.125, 25);

% plot_range_fft(ranges, params);
t = tiledlayout(1,2);
nexttile;
plot_range_doppler(ranges, velocities, params);
params.swerling = false;
nexttile;
plot_range_doppler(ranges, velocities, params);