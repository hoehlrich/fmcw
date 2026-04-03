function plot_cfar_2d(detected, z, T, params)
    t = 0:1/params.fs:params.Tc;
    N = length(t);

    f = (0:N)*(params.fs/N);
    R = (f*params.c*params.Tc)/(2*params.B);
    f_doppler = (-params.Nc/2:params.Nc/2-1)*(1/params.Tc)/params.Nc;
    v_axis = f_doppler*params.c/(2*params.fc);

    [r_idx, v_idx] = find(detected);

    %figure;
    subplot(1, 2, 2);
    imagesc(v_axis, R(1:N/2), 20*log10(abs(T)));
    colormap(jet);
    colorbar;
    axis xy;
    xlabel('Velocity (m/s)');
    ylabel('Range (m)');
    title('2D CFAR');

    hold on;
    scatter(v_axis(v_idx), R(r_idx), 'w+', 'LineWidth', 1.5);
end
