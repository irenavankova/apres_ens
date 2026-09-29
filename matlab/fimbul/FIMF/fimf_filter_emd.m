function [iNaN, iNaN_bed, iNaN_int] = fimf_filter_emd(~, ~, t_bed, t_int, y_bed_raw, y_bed_std, y_int_raw, y_int_std, ymax, window_len, ref_len, quant_grid, plot_results)
%FIMF_FILTER_EMD Apply 1D Wasserstein/Earth Mover's Distance filtering.
%   Inputs:
%       window_len - Length of local window (e.g., 81)
%       ref_len    - Length of reference window (e.g., 251)
%       quant_grid - Quantile evaluation grid (e.g., linspace(0.05, 0.95, 20))

    half_w = floor(window_len / 2);
    half_ref = floor(ref_len / 2);
    N_bed = length(y_bed_raw);
    N_int = length(y_int_raw);

    % Initialize outputs
    y_bed = y_bed_raw;
    y_int = y_int_raw;

    % Process bed data
    if N_bed >= window_len + 2 * ref_len
        emd_dist = zeros(N_bed, 1);
        for i = (half_ref + 1):(N_bed - half_ref)
            local_win = y_bed_raw(i - half_w : i + half_w);
            ref_win = y_bed_raw([i - half_ref : i - half_w, i + half_w : i + half_ref]);
            l_quantiles = quantile(local_win, quant_grid);
            r_quantiles = quantile(ref_win, quant_grid);
            emd_dist(i) = mean(abs(l_quantiles - r_quantiles));
        end
        thresh_emd = median(emd_dist, 'omitnan') + 4.5 * mad(emd_dist, 1);
        mask_emd = emd_dist > thresh_emd;
        y_bed(mask_emd) = NaN;
    end

    % Process internal data
    if N_int >= window_len + 2 * ref_len
        emd_dist = zeros(N_int, 1);
        for i = (half_ref + 1):(N_int - half_ref)
            local_win = y_int_raw(i - half_w : i + half_w);
            ref_win = y_int_raw([i - half_ref : i - half_w, i + half_w : i + half_ref]);
            l_quantiles = quantile(local_win, quant_grid);
            r_quantiles = quantile(ref_win, quant_grid);
            emd_dist(i) = mean(abs(l_quantiles - r_quantiles));
        end
        thresh_emd = median(emd_dist, 'omitnan') + 4.5 * mad(emd_dist, 1);
        mask_emd = emd_dist > thresh_emd;
        y_int(mask_emd) = NaN;
    end

    % Compute individual and combined NaN masks
    iNaN_bed = isnan(y_bed);
    iNaN_int = isnan(y_int);
    iNaN = iNaN_bed | iNaN_int;

    % Plot results if enabled
    if plot_results
        Lwide = 1;

        % Reconstruct combo for plotting
        y_bed_combo_plot = y_bed_raw;
        y_bed_combo_plot(iNaN) = NaN;
        y_int_combo_plot = y_int_raw;
        y_int_combo_plot(iNaN) = NaN;

        % --- FIGURE: Raw and filtered timeseries ---
        figure;
        subplot(2, 1, 1); h = gca;
        plot(t_bed(1:end-1), y_bed_raw, 'k-', 'LineWidth', Lwide); hold on;
        plot(t_bed(1:end-1), y_bed, 'r-', 'LineWidth', Lwide);
        plot(t_bed(1:end-1), y_bed_combo_plot, 'c-', 'LineWidth', Lwide);
        grid on;
        ylabel('Bed');
        datetick('x', 'mm/yy', 'keeplimits');
        title(['EMD Filter, window = ' num2str(window_len) ', ref = ' num2str(ref_len)], 'fontweight', 'normal');

        subplot(2, 1, 2); h = [h; gca];
        plot(t_int(1:end-1), y_int_raw, 'k-', 'LineWidth', Lwide); hold on;
        plot(t_int(1:end-1), y_int, 'r-', 'LineWidth', Lwide);
        plot(t_int(1:end-1), y_int_combo_plot, 'c-', 'LineWidth', Lwide);
        linkaxes(h, 'x');
        datetick('x', 'mm/yy', 'keeplimits');
        grid on;
        ylabel('Internal');

        % --- FIGURE: Uncertainty bounds ---
        figure;
        subplot(2, 1, 1); h = gca;
        plot(t_bed(1:end-1), y_bed_raw - y_bed_std, 'c-', 'LineWidth', Lwide); hold on;
        plot(t_bed(1:end-1), y_bed_raw + y_bed_std, 'c-', 'LineWidth', Lwide);
        plot(t_bed(1:end-1), y_bed, 'r-', 'LineWidth', Lwide);
        %plot(t_bed(1:end-1), y_bed_combo_plot, 'k-', 'LineWidth', Lwide);
        grid on;
        ylabel('Bed');
        datetick('x', 'mm/yy', 'keeplimits');
        ylim([-ymax ymax]);
        title(['EMD Filter, window = ' num2str(window_len) ', ref = ' num2str(ref_len)], 'fontweight', 'normal');

        subplot(2, 1, 2); h = [h; gca];
        plot(t_int(1:end-1), y_int_raw - y_int_std, 'c-', 'LineWidth', Lwide); hold on;
        plot(t_int(1:end-1), y_int_raw + y_int_std, 'c-', 'LineWidth', Lwide);
        plot(t_int(1:end-1), y_int, 'r-', 'LineWidth', Lwide);
        %plot(t_int(1:end-1), y_int_combo_plot, 'k-', 'LineWidth', Lwide);
        linkaxes(h, 'x');
        datetick('x', 'mm/yy', 'keeplimits');
        grid on;
        ylabel('Internal');
        ylim([-ymax ymax]);
    end
end