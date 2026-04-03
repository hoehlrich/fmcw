function plot_cfar_1d(detected, z, T, params)
    t = 0:1/params.fs:params.Tc;
    f = (0:length(t)/2-1)*(params.fs/length(t));
    R = (f*params.c*params.Tc)/(2*params.B);
    figure
    plot(R, 10*log10(T));
    hold on;
    z_dB = 10*log10(z);
    plot(R, z_dB);
    scatter(R(detected), z_dB(detected), 'rv', 'filled');
    legend('Threshold', 'Power');
    xlabel('Range');
    ylabel('dB');
    ylim([mean(z_dB(isfinite(z_dB))), max(z_dB(isfinite(z_dB))) + 5]);
    title('CFAR');
    grid on;
end
