function plot_range_doppler(ranges, velocities, params)
    c = 3e8;    % Speed of light
    d_step = velocities*params.Tc;

    t = 0:1/params.fs:params.Tc;
    N = length(t); % Num samples per chirp = num range fft bins

    % Populate range_profiles from ranges and velocities
    range_profiles = zeros(params.Nc, N);

    % Handle swerling if applicable
    if params.swerling
        A = raylrnd(1, length(ranges), params.Nc) / sqrt(pi/2);
    else
        A = ones(length(ranges), params.Nc);
    end

    for i = 1:params.Nc
        IF = zeros(size(t));
        for j = 1:length(ranges)
            tau = 2*ranges(j)/c;
            f_if = (params.B/params.Tc)*tau;
            phi = 2*pi*(2*ranges(j)/c)*params.fc;
            IF = IF + A(j,i)*exp(1j*(2*pi*f_if*t + phi));
        end
        IF = IF + params.noise*(randn(size(IF)) + 1j*randn(size(IF)));
        range_profiles(i, :) = fft(IF);
        ranges = ranges + d_step;
    end

    % Calculate range doppler matrix
    range_doppler = fftshift(fft(range_profiles, [], 1), 1);

    half_N = floor(N/2);
    f = (0:N)*(params.fs/N);
    R = (f*c*params.Tc)/(2*params.B);
    f_doppler = (-params.Nc/2:params.Nc/2-1)*(1/params.Tc)/params.Nc;
    v_axis = f_doppler*c/(2*params.fc);


    % figure;
    imagesc(v_axis, R(1:half_N), abs(range_doppler(:, 1:half_N))');
    colormap(jet);
    clim([0 0.1*max(abs(range_doppler(:)))]);
    axis xy;
    xlabel('Velocity (m/s)');
    ylabel('Range (m)');
   title('Range-Doppler Map');
    colorbar;
end
