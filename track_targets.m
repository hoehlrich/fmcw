% track_targets: runs Kalman tracking on objects for params.tracking.N_frames frames
%   history.true_positions    N_frames x N_targets
%   history.true_velocities   N_frames x N_targets
%   history.detections        cell array, one set of (range,vel) pairs per frame
%   history.track_states      N_frames x N_tracks x 2
%   history.track_covariances N_frames x N_tracks x 2x2
function history = track_targets(ranges, velocities, params)
    tracking = params.tracking;
    N_tracks = length(ranges);

    % Conversions
    t = 0:1/params.fs:params.Tc;
    N = length(t);
    f = (0:N)*(params.fs/N);
    R = (f*params.c*params.Tc)/(2*params.B);
    f_doppler = (-params.Nc/2:params.Nc/2-1)*(1/params.Tc)/params.Nc;
    v_axis = f_doppler*params.c/(2*params.fc);
    delta_t = params.Nc*params.Tc;

    P0 = tracking.R;
    for i = 1:N_tracks
        tracks(i).x = [ranges(i); velocities(i)];  % initial state from true positions
        tracks(i).P = P0;
    end
    for frame = 1:tracking.N_frames
        % generate range-Doppler map
        range_doppler_map = range_doppler(ranges, velocities, params);

        % run CFAR to get detections
        N = floor(params.Tc/(1/params.fs));
        [detected, ~, ~] = cfar_2d(range_doppler_map(:, 1:N/2)', params);
        detections = get_detections(detected, R, v_axis);
        history.detections{frame} = detections;

        for i = 1:N_tracks
            % update history
            history.true_positions(frame, i) = ranges(i);
            history.true_velocities(frame, i) = velocities(i);
            history.track_states(frame, i, :) = tracks(i).x;
            history.track_covariances(frame, i, :, :) = tracks(i).P;


            % run Kalman predict
            x_pred = tracking.F*tracks(i).x;
            P_pred = tracking.F*tracks(i).P*tracking.F' + tracking.Q;

            % associate detections to tracks (Mahalanobis gating)
            in_gate = false(size(detections, 1), 1);
            S = tracking.H*P_pred*tracking.H' + tracking.R;
            for k = 1:size(detections, 1)
                d = detections(k,:)' - x_pred;
                if d'*(S\d) <= tracking.gamma
                    in_gate(k) = true;
                end
            end
            z = mean(detections(in_gate, :), 1)';  % centroid [range; velocity]

            if isempty(detections) || sum(in_gate) == 0
                % no associated detections just use prediction
                tracks(i).x = x_pred;
                tracks(i).P = P_pred;
            else
                % there are associated detections do Kalman update
                K = P_pred*tracking.H'*(tracking.H*P_pred*tracking.H'+tracking.R)^(-1);
                x_est = x_pred + K*(z - tracking.H*x_pred);
                P_est = (eye(2) - K*tracking.H)*P_pred*(eye(2) - K*tracking.H)' + K*tracking.R*K';
                tracks(i).x = x_est;
                tracks(i).P = P_est;
            end

            % update target position
            ranges(i) = ranges(i) + velocities(i)*delta_t;
        end
    end
end

function detections = get_detections(detected, R, v_axis)
    [r_idx, v_idx] = find(detected);
    detections = [R(r_idx)', v_axis(v_idx)'];
end
