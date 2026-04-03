function plot_range_doppler(range_doppler_map, params)
    t = 0:1/params.fs:params.Tc;
    N = length(t); % Num samples per chirp = num range fft bins

    half_N = floor(N/2);
    f = (0:N)*(params.fs/N);
    R = (f*params.c*params.Tc)/(2*params.B);
    f_doppler = (-params.Nc/2:params.Nc/2-1)*(1/params.Tc)/params.Nc;
    v_axis = f_doppler*params.c/(2*params.fc);

    power = 20*log10(abs(range_doppler_map(:, 1:half_N)'));

    % figure;
    imagesc(v_axis, R(1:half_N), power);
    clim([mean(power(:)), max(power(:))]);
    colormap(jet);
    colorbar;
    axis xy;
    xlabel('Velocity (m/s)');
    ylabel('Range (m)');
    title('Range-Doppler Map');

    % [detected, ~, ~] = cfar_2d(range_doppler_map(:, 1:half_N)', params);
    % [r_idx, v_idx] = find(detected);
    % hold on;
    % scatter(v_axis(v_idx), R(r_idx), 'w+', 'LineWidth', 1);
    % hold off;

end
