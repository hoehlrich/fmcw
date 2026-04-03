% cfar_1d: return array 'detected' full of indices of cells that were
%          detetected for given params 'cfar' along with power and
%          threshold vectors
function [detected, z, T] = cfar_1d(cells, params)
    cfar = params.cfar_1d;
    N = length(cells);
    alpha = cfar.N_R*(cfar.P_FA^(-1/cfar.N_R) - 1); % CA-CFAR constant
    assert(N > cfar.N_R + cfar.N_G + 1, "Invalid parameters");
    detected = [];
    T = zeros(N);
    z = zeros(N);
    for i = 1:N
        left_overlap = max(-(i - cfar.N_G/2 - 1 - cfar.N_R), 0);
        right_overlap = max((i + 1 + cfar.N_G/2 + cfar.N_R/2) - N, 0);
        lag_head = i - cfar.N_G/2 - 1;  % Index of the lagging cell closest to CUT
        lead_head = i + cfar.N_G/2 + 1; % Index of the leading cell closest to CUT

        lag_length = cfar.N_R/2 - left_overlap + right_overlap;
        lead_length = cfar.N_R/2 - right_overlap + left_overlap;

        y = [cells(lag_head-lag_length:lag_head), cells(lead_head:lead_head+lead_length)];
        sigma = sum(y.^2)/cfar.N_R; % Mean interference power
        T(i) = alpha*sigma;
        z(i) = cells(i)^2;
        if z(i) > T(i)
            detected = [detected, i];
        end
    end
end
