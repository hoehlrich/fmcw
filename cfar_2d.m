% cfar_2d: return logical matrix 'detected' of detections for given params
%          along with power and threshold matrices
function [detected, z, T] = cfar_2d(map, params)
    cfar = params.cfar_2d;
    [N_range, N_vel] = size(map);

    % Total reference cells (full outer window minus guard+CUT window)
    N_ref = (2*cfar.N_R + 2*cfar.N_G + 1) * (2*cfar.N_CR + 2*cfar.N_CG + 1) ...
          - (2*cfar.N_G + 1) * (2*cfar.N_CG + 1);

    assert(N_range > 2*(cfar.N_R + cfar.N_G) + 1, "Invalid range parameters");
    assert(N_vel   > 2*(cfar.N_CR + cfar.N_CG) + 1, "Invalid cross-range parameters");

    z = abs(map).^2; % Square law detection
    T = zeros(N_range, N_vel);
    detected = false(N_range, N_vel);

    for ri = 1:N_range
        for vi = 1:N_vel
            % Outer window bounds (reference + guard + CUT)
            r_out_min = ri - cfar.N_R - cfar.N_G;
            r_out_max = ri + cfar.N_R + cfar.N_G;
            v_out_min = vi - cfar.N_CR - cfar.N_CG;
            v_out_max = vi + cfar.N_CR + cfar.N_CG;

            % Inner window bounds (guard cells + CUT)
            r_in_min = ri - cfar.N_G;
            r_in_max = ri + cfar.N_G;
            v_in_min = vi - cfar.N_CG;
            v_in_max = vi + cfar.N_CG;

            % Clamp outer window to map boundaries
            r_out_min_c = max(r_out_min, 1);
            r_out_max_c = min(r_out_max, N_range);
            v_out_min_c = max(v_out_min, 1);
            v_out_max_c = min(v_out_max, N_vel);

            % Extract outer window
            outer = z(r_out_min_c:r_out_max_c, v_out_min_c:v_out_max_c);

            % Map guard region into clamped outer window coordinates
            g_r_min = r_in_min - r_out_min_c + 1;
            g_r_max = r_in_max - r_out_min_c + 1;
            g_v_min = v_in_min - v_out_min_c + 1;
            g_v_max = v_in_max - v_out_min_c + 1;

            % Clamp guard mask to outer window size
            g_r_min = max(g_r_min, 1);
            g_r_max = min(g_r_max, size(outer, 1));
            g_v_min = max(g_v_min, 1);
            g_v_max = min(g_v_max, size(outer, 2));

            % Reference mask: outer window minus guard+CUT region
            mask = true(size(outer));
            mask(g_r_min:g_r_max, g_v_min:g_v_max) = false;

            ref_cells = outer(mask);
            n_ref_actual = length(ref_cells);

            if n_ref_actual == 0
                continue;
            end

            % Recompute alpha locally to handle edge cells correctly
            alpha_local = n_ref_actual * (cfar.P_FA^(-1/n_ref_actual) - 1);

            sigma = sum(ref_cells) / n_ref_actual;
            T(ri, vi) = alpha_local * sigma;

            if z(ri, vi) > T(ri, vi)
                detected(ri, vi) = true;
            end
        end
    end
end
