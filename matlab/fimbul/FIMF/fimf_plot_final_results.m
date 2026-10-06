function fimf_plot_final_results(t_bed, y_bed_merged_yr, y_bed_ci_yr, vsr_lin, vsr_lin_se, ...
                                 vvel_t, vvel_vsr_lin, vvel_vsr_lin_se, vvel_v_bed, vvel_v_bed_se)
% FIMF_PLOT_FINAL_RESULTS Generates final summary figures for VSR and Basal Melt Rate.

    % Ensure column vectors
    t_bed = t_bed(:);
    vvel_t = vvel_t(:);
    
    %% Figure 1: Sliding Window Summary
    figure('Name', 'Sliding Window Summary', 'Position', [100, 100, 800, 600]);

    % --- VSR Interpolation with END-POINT Extrapolation ---
    valid_interp_idx = ~isnan(vvel_vsr_lin);
    t_valid = vvel_t(valid_interp_idx);
    vsr_valid = vvel_vsr_lin(valid_interp_idx);
    vsr_se_valid = vvel_vsr_lin_se(valid_interp_idx);

    vsr_interp = interp1(t_valid, vsr_valid, t_bed, 'linear');
    vsr_se_interp = interp1(t_valid, vsr_se_valid, t_bed, 'linear');

    % Hold values constant outside the bounds
    vsr_interp(t_bed <= t_valid(1)) = vsr_valid(1);
    vsr_interp(t_bed >= t_valid(end)) = vsr_valid(end);
    vsr_se_interp(t_bed <= t_valid(1)) = vsr_se_valid(1);
    vsr_se_interp(t_bed >= t_valid(end)) = vsr_se_valid(end);

    % --- Weighted Linear Fit to VSR Over Time ---
    t_mean = mean(vvel_t, 'omitnan');
    t_fit = vvel_t - t_mean;
    vsr_fit_data = vvel_vsr_lin;
    w_fit = 1 ./ (vvel_vsr_lin_se.^2); 

    valid_fit_idx = ~isnan(vsr_fit_data) & ~isnan(w_fit) & ~isinf(w_fit);
    t_fit_valid = t_fit(valid_fit_idx);
    vsr_fit_valid = vsr_fit_data(valid_fit_idx);
    w_fit_valid = w_fit(valid_fit_idx);

    A_fit = [t_fit_valid, ones(length(t_fit_valid), 1)];
    [p_vsr, ~] = lscov(A_fit, vsr_fit_valid, w_fit_valid);

    t_bed_centered = t_bed - t_mean;
    vsr_trend = p_vsr(1) * t_bed_centered + p_vsr(2);
    vsr_trend_se = mean(vvel_vsr_lin_se, 'omitnan'); 

    % Subplot 1
    subplot(2,1,1); h_fig1 = gca; hold on;
    errorbar(t_bed, y_bed_merged_yr, y_bed_ci_yr, 'Color', [0.8 0.8 0.8], 'LineStyle', 'none', ...
        'CapSize', 0, 'DisplayName', '95% CI (Raw)');
    plot(t_bed, y_bed_merged_yr, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, ...
        'DisplayName', 'y\_bed\_merged'); 
    errorbar(vvel_t, vvel_v_bed, vvel_v_bed_se, 'r.-', 'MarkerSize', 12, 'LineWidth', 1, ...
        'DisplayName', 'Bed Rate of Change (Sliding)');
    datetick('x', 'mm/yy', 'keeplimits');
    ylabel('Velocity (m/a)'); title('Basal Reflector Velocity Comparison');
    legend('Location', 'best'); grid on;

    % Subplot 2
    subplot(2,1,2); h_fig1 = [h_fig1; gca]; hold on;
    
    ci_lower_vsr = vsr_lin - vsr_lin_se;
    ci_upper_vsr = vsr_lin + vsr_lin_se;
    x_fill_vsr = [t_bed(:); flipud(t_bed(:))];
    y_fill_vsr = [t_bed(:)*0 + ci_lower_vsr; t_bed(:)*0 + ci_upper_vsr];
    fill(x_fill_vsr, y_fill_vsr, 0.5*[1 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2, ...
        'DisplayName', 'Overall Mean \pm SE');

    ci_lower_trend = vsr_trend - vsr_trend_se;
    ci_upper_trend = vsr_trend + vsr_trend_se;
    y_fill_trend = [ci_lower_trend(:); flipud(ci_upper_trend(:))];
    fill(x_fill_vsr, y_fill_trend, [0.4 0.8 0.4], 'EdgeColor', 'none', 'FaceAlpha', 0.3, ...
        'DisplayName', 'Linear Fit \pm SE');

    plot(t_bed, t_bed*0 + vsr_lin, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Overall Mean VSR');
    plot(t_bed, vsr_trend, 'g--', 'LineWidth', 2, 'DisplayName', 'Linear Fit VSR');
    plot(t_bed, vsr_interp, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Extrapolated VSR (End-point)');
    errorbar(vvel_t, vvel_vsr_lin, vvel_vsr_lin_se, 'b.-', 'MarkerSize', 12, 'LineWidth', 1, ...
        'DisplayName', 'Linear VSR (Sliding)');
    datetick('x', 'mm/yy', 'keeplimits');
    ylabel('VSR (a^{-1})'); title('Linear Vertical Strain Rate over Time');
    legend('Location', 'best'); grid on;
    linkaxes(h_fig1, 'x');


    %% Figure 2: Basal Melt Rate Figure (4 Panels)
    figure('Name', 'Basal Melt Rate Summary', 'Position', [150, 150, 800, 1000]);

    bmr_const = -1 * (y_bed_merged_yr - vsr_lin);
    bmr_const_err = sqrt(y_bed_ci_yr.^2 + vsr_lin_se^2);

    bmr_trend = -1 * (y_bed_merged_yr - vsr_trend);
    bmr_trend_err = sqrt(y_bed_ci_yr.^2 + vsr_trend_se.^2);

    bmr_interp = -1 * (y_bed_merged_yr - vsr_interp);
    bmr_interp_err = sqrt(y_bed_ci_yr.^2 + vsr_se_interp.^2);

    bmr_slide = -1 * (vvel_v_bed - vvel_vsr_lin);
    bmr_slide_err = sqrt(vvel_v_bed_se.^2 + vvel_vsr_lin_se.^2);

    % Panel 1
    subplot(4,1,1); h_fig2 = gca; hold on;
    errorbar(t_bed, bmr_const, bmr_const_err, 'Color', [0.8 0.8 0.8], 'LineStyle', 'none', ...
        'CapSize', 0, 'DisplayName', 'Geometric Uncertainty');
    plot(t_bed, bmr_const, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, 'DisplayName', 'BMR (Constant VSR)');
    datetick('x', 'mm/yy', 'keeplimits'); ylabel('BMR (m/a)');
    title('Basal Melt Rate (High-Res Bed, Constant Overall VSR)'); legend('Location', 'best'); grid on;

    % Panel 2
    subplot(4,1,2); h_fig2 = [h_fig2; gca]; hold on;
    errorbar(t_bed, bmr_trend, bmr_trend_err, 'Color', [0.8 0.8 0.8], 'LineStyle', 'none', ...
        'CapSize', 0, 'DisplayName', 'Geometric Uncertainty');
    plot(t_bed, bmr_trend, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, 'DisplayName', 'BMR (Linear Fit VSR)');
    datetick('x', 'mm/yy', 'keeplimits'); ylabel('BMR (m/a)');
    title('Basal Melt Rate (High-Res Bed, Linear Fit VSR)'); legend('Location', 'best'); grid on;

    % Panel 3
    subplot(4,1,3); h_fig2 = [h_fig2; gca]; hold on;
    errorbar(t_bed, bmr_interp, bmr_interp_err, 'Color', [0.8 0.8 0.8], 'LineStyle', 'none', ...
        'CapSize', 0, 'DisplayName', 'Geometric Uncertainty');
    plot(t_bed, bmr_interp, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, 'DisplayName', 'BMR (End-point VSR)');
    datetick('x', 'mm/yy', 'keeplimits'); ylabel('BMR (m/a)');
    title('Basal Melt Rate (High-Res Bed, End-Point Extrapolated VSR)'); legend('Location', 'best'); grid on;

    % Panel 4
    subplot(4,1,4); h_fig2 = [h_fig2; gca]; hold on;
    errorbar(vvel_t, bmr_slide, bmr_slide_err, 'm.-', 'MarkerSize', 12, 'LineWidth', 1.5, ...
        'DisplayName', 'BMR (Sliding Window)');
    datetick('x', 'mm/yy', 'keeplimits'); ylabel('BMR (m/a)');
    title('Basal Melt Rate (Sliding Window Estimates)'); legend('Location', 'best'); grid on;
    linkaxes(h_fig2, 'xy');


    %% Figure 3: Overlay Basal Melt Rate Figure
    figure('Name', 'BMR Overlay Summary', 'Position', [200, 200, 1000, 500]); hold on;

    c_const  = [0.000, 0.447, 0.741]; 
    c_trend  = [0.466, 0.674, 0.188]; 
    c_interp = [0.850, 0.325, 0.098]; 
    x_fill = [t_bed(:); flipud(t_bed(:))];

    % 1. Constant VSR
    bmr_c_fill = bmr_const(:); err_c_fill = bmr_const_err(:);
    nan_idx_c = isnan(bmr_c_fill) | isnan(err_c_fill);
    err_c_fill(nan_idx_c) = 0; valid_c = ~nan_idx_c;
    if any(nan_idx_c) && any(valid_c), bmr_c_fill(nan_idx_c) = interp1(t_bed(valid_c), bmr_c_fill(valid_c), t_bed(nan_idx_c), 'linear', 'extrap'); end
    y_fill_const = [(bmr_c_fill - err_c_fill); flipud(bmr_c_fill + err_c_fill)];

    % 2. Linear Fit VSR
    bmr_t_fill = bmr_trend(:); err_t_fill = bmr_trend_err(:);
    nan_idx_t = isnan(bmr_t_fill) | isnan(err_t_fill);
    err_t_fill(nan_idx_t) = 0; valid_t = ~nan_idx_t;
    if any(nan_idx_t) && any(valid_t), bmr_t_fill(nan_idx_t) = interp1(t_bed(valid_t), bmr_t_fill(valid_t), t_bed(nan_idx_t), 'linear', 'extrap'); end
    y_fill_trend = [(bmr_t_fill - err_t_fill); flipud(bmr_t_fill + err_t_fill)];

    % 3. Interpolated/End-Point VSR
    bmr_i_fill = bmr_interp(:); err_i_fill = bmr_interp_err(:);
    nan_idx_i = isnan(bmr_i_fill) | isnan(err_i_fill);
    err_i_fill(nan_idx_i) = 0; valid_i = ~nan_idx_i;
    if any(nan_idx_i) && any(valid_i), bmr_i_fill(nan_idx_i) = interp1(t_bed(valid_i), bmr_i_fill(valid_i), t_bed(nan_idx_i), 'linear', 'extrap'); end
    y_fill_interp = [(bmr_i_fill - err_i_fill); flipud(bmr_i_fill + err_i_fill)];

    fill(x_fill, y_fill_const, c_const, 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off');
    fill(x_fill, y_fill_trend, c_trend, 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off');
    fill(x_fill, y_fill_interp, c_interp, 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off');

    plot(t_bed, bmr_const, '-', 'Color', c_const, 'LineWidth', 1.5, 'DisplayName', 'BMR (Constant VSR)');
    plot(t_bed, bmr_trend, '-', 'Color', c_trend, 'LineWidth', 1.5, 'DisplayName', 'BMR (Linear Fit VSR)');
    plot(t_bed, bmr_interp, '-', 'Color', c_interp, 'LineWidth', 1.5, 'DisplayName', 'BMR (End-point VSR)');

    datetick('x', 'mm/yy', 'keeplimits');
    ylabel('BMR (m/a)'); title('Basal Melt Rate Comparison Overlay');
    legend('Location', 'best'); grid on;
end