function plot_range_fft(ranges, params)
    c = 3e8;    % Speed of light
    t = 0:1/params.fs:params.Tc;

    % Generate IF signal from ranges
    IF = zeros(size(t));
    for i = 1:length(ranges)
        tau = 2*ranges(i)/c;
        f_if = (params.B/params.Tc)*tau;
        IF = IF + cos(2*pi*f_if*t);
    end
    IF = IF + 0.7*randn(size(IF));

    N = length(t);
    IF_fft = abs(fft(IF));
    IF_fft = IF_fft(1:N/2);
    f = (0:N/2-1)*(params.fs/N);
    R = (f*c*params.Tc)/(2*params.B);

    figure;
    plot(R, IF_fft);
    xlabel('Range');
    ylabel('Amplitude');
    title('Range FFT');
    grid on;
end
