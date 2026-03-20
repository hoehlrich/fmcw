function [new_ranges, new_velocities] = spread_targets(ranges, velocities, spread_r, spread_v, n)
    new_ranges = [];
    new_velocities = [];
    for i = 1:length(ranges)
        r_points = [ranges(i), ranges(i) + (rand(1,n)-0.5)*2*spread_r];
        v_points = [velocities(i), velocities(i) + (rand(1,n)-0.5)*2*spread_v];
        new_ranges = [new_ranges, r_points];
        new_velocities = [new_velocities, v_points];
    end
end
