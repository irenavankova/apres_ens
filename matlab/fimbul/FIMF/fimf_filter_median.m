function [iNaN, iNaN_bed, iNaN_int] = fimf_filter_median(dhdtmat_bed, dhdtmat_int, t_bed, t_int, y_bed_raw, y_bed_std, y_int_raw, y_int_std, ymax, plot_results)
%FIMF_FILTER_MEDIAN Apply median outlier detection.

% Initialize filtered outputs
y_bed = [];
y_int = [];

% Process bed data
if ~isempty(y_bed_raw)
    ii = isoutlier(y_bed_raw, 'median');
    ind = find(ii == 1);
    y_bed = y_bed_raw;
    y_bed(ind) = NaN;
    y_bed(ind+1) = NaN;
    y_bed(ind+2) = NaN;
    try y_bed(ind-2) = NaN; end
    try y_bed(ind-1) = NaN; end
end

% Process internal data
if ~isempty(y_int_raw)
    ii = isoutlier(y_int_raw, 'median');
    ind = find(ii == 1);
    y_int = y_int_raw;
    y_int(ind) = NaN;
    y_int(ind+1) = NaN;
    y_int(ind+2) = NaN;
    try y_int(ind-2) = NaN; end
    try y_int(ind-1) = NaN; end
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
    title('Median Filter - Raw and Filtered', 'fontweight', 'normal');

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
    title('Median Filter - Uncertainty Bounds', 'fontweight', 'normal');

    subplot(2, 1, 2); h = [h; gca];
    plot(t_int(1:end-1), y_int_raw - y_int_std, 'c-', 'LineWidth', Lwide); hold on;
    plot(t_int(1:end-1), y_int_raw + y_int_std, 'c-', 'LineWidth', Lwide);
    plot(t_int(1:end-1), y_int, 'r-', 'LineWidth', Lwide);
    plot(t_int(1:end-1), y_int_combo_plot, 'k-', 'LineWidth', Lwide);
    linkaxes(h, 'x');
    datetick('x', 'mm/yy', 'keeplimits');
    grid on;
    ylabel('Internal');
    ylim([-ymax ymax]);
end
end