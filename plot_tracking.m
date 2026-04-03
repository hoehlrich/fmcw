function plot_tracking(history, params)
    N_frames = params.tracking.N_frames;
    N_tracks = size(history.true_positions, 2);
    colors = lines(N_tracks);

    figure;
    subplot(2,1,1);
    hold on;
    for i = 1:N_tracks
        plot(1:N_frames, history.true_positions(:,i), '--', 'Color', colors(i,:), 'LineWidth', 1.5);
        plot(1:N_frames, squeeze(history.track_states(:,i,1)), '-', 'Color', colors(i,:), 'LineWidth', 1.5);
    end
    ylabel('Range (m)');
    xlabel('Frame');
    title('Range Tracking');
    legend('True', 'Estimated');
    grid on;

    subplot(2,1,2);
    hold on;
    for i = 1:N_tracks
        plot(1:N_frames, history.true_velocities(:,i), '--', 'Color', colors(i,:), 'LineWidth', 1.5);
        plot(1:N_frames, squeeze(history.track_states(:,i,2)), '-', 'Color', colors(i,:), 'LineWidth', 1.5);
    end
    ylabel('Velocity (m/s)');
    xlabel('Frame');
    title('Velocity Tracking');
    grid on;
end
