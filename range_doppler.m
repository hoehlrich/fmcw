function range_doppler_map = range_doppler(ranges, velocities, params)
    c = 3e8;    % Speed of light
    d_step = velocities*params.Tc;

    t = 0:1/params.fs:params.Tc;
    N = length(t); % Num samples per chirp = num range fft bins

    % Populate range_profiles from ranges and velocities
    range_profiles = zeros(params.Nc, N);

    % Handle swerling if applicable
    if params.swerling
        A = raylrnd(0.5, length(ranges), params.Nc) / sqrt(pi/2);
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
    range_doppler_map = fftshift(fft(range_profiles, [], 1), 1);
end
