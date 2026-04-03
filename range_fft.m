function IF_fft = range_fft(ranges, params)
    t = 0:1/params.fs:params.Tc;

    % Generate IF signal from ranges
    IF = zeros(size(t));
    for i = 1:length(ranges)
        tau = 2*ranges(i)/params.c;
        f_if = (params.B/params.Tc)*tau;
        IF = IF + cos(2*pi*f_if*t);
    end
    IF = IF + params.noise*randn(size(IF));

    N = length(t);
    IF_fft = abs(fft(IF));
    IF_fft = IF_fft(1:N/2);
end
