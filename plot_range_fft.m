function plot_range_fft(IF_fft, params)
    t = 0:1/params.fs:params.Tc;
    N = length(t);
    f = (0:N/2-1)*(params.fs/N);
    R = (f*params.c*params.Tc)/(2*params.B);

    figure;
    plot(R, IF_fft);
    if params.cfar
        hold on;
        [detected, ~, ~] = cfar_1d(IF_fft, params);
        scatter(R(detected), IF_fft(detected), 'rv', 'filled');
    end
    xlabel('Range');
    ylabel('Amplitude');
    title('Range FFT');
    grid on;
end
