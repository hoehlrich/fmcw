fs = 20e6;  % Sampling frequency
Tc = 40e-6; % Chirp period
B = 4e9;    % Chirp bandwidth
c = 3e8;    % Speed of light
t = 0:1/fs:Tc;

ranges = [5,6,12];
% Generate IF signal from ranges
IF = zeros(size(t));
for i = 1:length(ranges)
    tau = 2*ranges(i)/c;
    f_if = (B/Tc)*tau;
    IF = IF + cos(2*pi*f_if*t);
end
IF = IF + 0.1*randn(size(IF));

N = length(IF);
IF_fft = abs(fft(IF));
IF_fft = IF_fft(1:N/2);
f = (0:N/2-1)*(fs/N);
R = (f*c*Tc)/(2*B);

figure;
plot(R, IF_fft);
xlabel('Range');
ylabel('Amplitude');
title('Range FFT');
grid on;
