function [iNaN, iNaN_bed, iNaN_int] = fimf_filter_mad(~, ~, t_bed, t_int, y_bed_raw, y_bed_std, y_int_raw, y_int_std, ymax, window_len, plot_results)
%FIMF_FILTER_MAD Apply Local Moving Dispersion (Rolling MAD) filtering.
%   Inputs:
%       window_len - Window length for movmad (e.g., 81)

% Initialize outputs
y_bed = y_bed_raw;
y_int = y_int_raw;

% Process bed data
if length(y_bed_raw) >= window_len
    mov_mad = movmad(y_bed_raw, window_len);
    thresh_mad = median(mov_mad, 'omitnan') + 5 * mad(mov_mad, 1);
    mask_mad = mov_mad > thresh_mad;
    y_bed(mask_mad) = NaN;
end

% Process internal data
if length(y_int_raw) >= window_len
    mov_mad = movmad(y_int_raw, window_len);
    thresh_mad = median(mov_mad, 'omitnan') + 5 * mad(mov_mad, 1);
    mask_mad = mov_mad > thresh_mad;
    y_int(mask_mad) = NaN;
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
    title(['MAD Filter, window = ' num2str(window_len)], 'fontweight', 'normal');

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
    plot(t_bed(1:end-1), y_bed_combo_plot, 'k-', 'LineWidth', Lwide);
    grid on;
    ylabel('Bed');
    datetick('x', 'mm/yy', 'keeplimits');
    ylim([-ymax ymax]);
    title(['MAD Filter, window = ' num2str(window_len)], 'fontweight', 'normal');

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