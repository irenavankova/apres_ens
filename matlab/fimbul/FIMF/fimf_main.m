%% Main Script: Sequential Execution of Filter Jumps and VSR Estimation
% Load data once, run filter jumps in 'final' mode to get iNaN_array,
% then pass to VSR estimation for both the whole timeseries and a sliding window.

%% Parameters
% Site and frequency settings
sitename = 'FIMFY25_B';
f0 = 200;
fend = 360;
db = 10;
bw = 40;
i1t = 3;

d2y = 365.25; 

% Range and filtering
dhRange_low = 20;
dhRange_high = 80;
opt_bed_source = 'tp'; % 'tp' (loads bed_tpq) or 'xcor' (loads bed_xcor)

% Filter parameters
Nbad = 6;
Nbaddh = 3;
ymax = 0.06;
plot_all_filters = false;
plot_results = true;

% VSR parameters
opt_combine_method = 'median'; 
opt_remove_outliers_vel = true; 
opt_vel_method = 'deriv_simple'; 
opt_vel_method_deriv_bed = 'mean'; 
opt_vel_method_deriv_int = 'median'; 
ylim_val = [0 170];
ytick_dist = 20;
drnf_vsr = [0 40; 124 4000];
drnf_artefacts = [];
max_nan_threshold = 0.5; % More than 20% bad -> do NaN

% Sliding Window parameters (in days)
%dt_window = 30; 
%dt_over = 5;
dt_window = 30; 
dt_over = 1;

%% Load all raw data once
[bed_data, int_data, t_bed, t_int, f1] = fimf_load_matfiles(...
    sitename, f0, fend, db, bw, i1t, dhRange_low, dhRange_high, opt_bed_source);

%% Run filter jumps in 'final' mode to get iNaN_array and uncertainties[cite: 20, 21]
[iNaN_array, y_bed_merged, y_int_merged, y_bed_std, y_bed_ci, y_int_std, y_int_ci] = fimf_filter_jumps(...
    bed_data, int_data, t_bed, t_int, f1, ...
    'sitename', sitename, ...
    'bw', bw, ...
    'Nbad', Nbad, ...
    'Nbaddh', Nbaddh, ...
    'i1t', i1t, ...
    'dhRange_low', dhRange_low, ...
    'dhRange_high', dhRange_high, ...
    'ymax', ymax, ...
    'plot_all_filters', plot_all_filters, ...
    'plot_results', plot_results, ...
    'mode', 'final', ...
    'opt_bed_source', opt_bed_source);

%% 1. Run VSR estimation over the WHOLE timeseries
plot_vsr_results_whole = true; 

disp('--- Calculating VSR over the WHOLE time series ---');
[vsr_lin, vsr_lin_se, vsr_quad, vsr_quad_se, v_bed, v_bed_se, z_bed_comb] = fimf_vsr_multifreq_estimation(...
    sitename, bed_data, int_data, t_bed, t_int, f1, d2y, ...
    'iNaN_array', iNaN_array, ...
    'opt_combine_method', opt_combine_method, ...
    'opt_remove_outliers_vel', opt_remove_outliers_vel, ...
    'opt_vel_method', opt_vel_method, ...
    'opt_vel_method_deriv_bed', opt_vel_method_deriv_bed, ...
    'opt_vel_method_deriv_int', opt_vel_method_deriv_int, ...
    'opt_bed_source', opt_bed_source, ...
    'ylim_val', ylim_val, ...
    'ytick_dist', ytick_dist, ...
    'drnf_vsr', drnf_vsr, ...
    'drnf_artefacts', drnf_artefacts, ...
    'max_nan_threshold', max_nan_threshold, ... 
    'plot_results', plot_vsr_results_whole);

% Display overall mean results
fprintf('--- VSR Estimation Results at Bed Depth (z = %.2f m) ---\n', z_bed_comb);
fprintf('Bed Rate of Change: %.4f +/- %.4f m/a\n', v_bed, v_bed_se);
fprintf('Linear VSR thinning rate:         %.6f +/- %.6f m/a\n', vsr_lin, vsr_lin_se);
fprintf('Quadratic VSR thinning rate:      %.6f +/- %.6f m/a\n', vsr_quad, vsr_quad_se);

%% 2. Run VSR estimation over SLIDING WINDOW
plot_vsr_results_sliding = false; 

disp('--- Calculating VSR over SLIDING WINDOW ---');

t_start_loop = t_int(1);
t_end_loop = t_int(end);

t1_win = t_start_loop;
t2_win = t1_win + dt_window;

ctr = 0;
vvel = struct();
vvel.t = []; vvel.v_bed = []; vvel.v_bed_se = []; 
vvel.vsr_lin = []; vvel.vsr_lin_se = [];
vvel.vsr_quad = []; vvel.vsr_quad_se = [];

% Prepare axes for sliding window animation (reused sequentially)
if plot_vsr_results_sliding
    [hax_anim, ~] = ap_sub_sub(2,1,3,10,0.25,0);
else
    hax_anim = [];
end

while t2_win <= t_end_loop
    ctr = ctr + 1;
    t_mid = t1_win + (t2_win - t1_win) / 2;
    
    [v_l, v_l_se, v_q, v_q_se, v_b, v_b_se, ~] = fimf_vsr_multifreq_estimation(...
        sitename, bed_data, int_data, t_bed, t_int, f1, d2y, ...
        'iNaN_array', iNaN_array, ...
        'opt_combine_method', opt_combine_method, ...
        'opt_remove_outliers_vel', opt_remove_outliers_vel, ...
        'opt_vel_method', opt_vel_method, ...
        'opt_vel_method_deriv_bed', opt_vel_method_deriv_bed, ...
        'opt_vel_method_deriv_int', opt_vel_method_deriv_int, ...
        'opt_bed_source', opt_bed_source, ...
        'ylim_val', ylim_val, ...
        'ytick_dist', ytick_dist, ...
        'drnf_vsr', drnf_vsr, ...
        'drnf_artefacts', drnf_artefacts, ...
        'max_nan_threshold', max_nan_threshold, ... 
        'plot_results', plot_vsr_results_sliding, ...
        't_start', t1_win, ...
        't_end', t2_win, ...
        'hax', hax_anim); % Pass the persistent axes to overwrite
        
    vvel.t(ctr) = t_mid;
    vvel.v_bed(ctr) = v_b;
    vvel.v_bed_se(ctr) = v_b_se;
    vvel.vsr_lin(ctr) = v_l;
    vvel.vsr_lin_se(ctr) = v_l_se;
    vvel.vsr_quad(ctr) = v_q;
    vvel.vsr_quad_se(ctr) = v_q_se;
    
    if plot_vsr_results_sliding
        drawnow;
    end
    
    t1_win = t1_win + dt_over;
    t2_win = t1_win + dt_window;
end

%% 3. Final Summary Figure
figure('Name', 'Sliding Window Summary', 'Position', [100, 100, 800, 600]);

% Convert arrays to m/a units for consistency
y_bed_merged_yr = y_bed_merged * d2y;
y_bed_ci_yr = y_bed_ci * d2y;

% Setup plotting limits (ensure column vector)
t_plot = reshape(t_bed(1:length(y_bed_merged_yr)), [], 1);

% --- VSR Interpolation with END-POINT Extrapolation ---
% Only interpolate using valid (non-NaN) sliding window points
valid_interp_idx = ~isnan(vvel.vsr_lin);
t_valid = vvel.t(valid_interp_idx);
vsr_valid = vvel.vsr_lin(valid_interp_idx);
vsr_se_valid = vvel.vsr_lin_se(valid_interp_idx);

vsr_interp = interp1(t_valid, vsr_valid, t_plot, 'linear');
vsr_se_interp = interp1(t_valid, vsr_se_valid, t_plot, 'linear');

% Hold values constant outside the bounds of the valid sliding window
vsr_interp(t_plot <= t_valid(1)) = vsr_valid(1);
vsr_interp(t_plot >= t_valid(end)) = vsr_valid(end);
vsr_se_interp(t_plot <= t_valid(1)) = vsr_se_valid(1);
vsr_se_interp(t_plot >= t_valid(end)) = vsr_se_valid(end);

% --- Weighted Linear Fit to VSR Over Time ---
% Center time to improve numerical stability in fitting
t_mean = mean(vvel.t, 'omitnan');
t_fit = vvel.t(:) - t_mean;
vsr_fit_data = vvel.vsr_lin(:);
w_fit = 1 ./ (vvel.vsr_lin_se(:).^2); % Weights based on sliding window SE

% Filter out NaNs and Infs before passing to lscov
valid_fit_idx = ~isnan(vsr_fit_data) & ~isnan(w_fit) & ~isinf(w_fit);
t_fit_valid = t_fit(valid_fit_idx);
vsr_fit_valid = vsr_fit_data(valid_fit_idx);
w_fit_valid = w_fit(valid_fit_idx);

A_fit = [t_fit_valid, ones(length(t_fit_valid), 1)];

% Calculate linear fit and coefficient errors
[p_vsr, se_vsr] = lscov(A_fit, vsr_fit_valid, w_fit_valid);

% Evaluate trend line and calculate geometric error band over timeline
t_plot_centered = t_plot - t_mean;
vsr_trend = p_vsr(1) * t_plot_centered + p_vsr(2);
vsr_trend_se = sqrt((se_vsr(1) * t_plot_centered).^2 + se_vsr(2)^2);


% Subplot 1: Bed Rate of Change vs y_bed_merged
subplot(2,1,1); h_fig1 = gca;
hold on;
% Plot raw CI as light gray vertical error bars without caps
errorbar(t_plot, y_bed_merged_yr, y_bed_ci_yr, 'Color', [0.8 0.8 0.8], 'LineStyle', 'none', ...
    'CapSize', 0, 'DisplayName', '95% CI (Raw)');
% Plot merged timeline and sliding window estimates on top
plot(t_plot, y_bed_merged_yr, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, ...
    'DisplayName', 'y\_bed\_merged'); 
errorbar(vvel.t, vvel.v_bed, vvel.v_bed_se, 'r.-', 'MarkerSize', 12, 'LineWidth', 1, ...
    'DisplayName', 'Bed Rate of Change (Sliding)');

datetick('x', 'mm/yy', 'keeplimits');
ylabel('Velocity (m/a)');
title('Basal Reflector Velocity Comparison');
legend('Location', 'best');
grid on;


% Subplot 2: Linear VSR
subplot(2,1,2); h_fig1 = [h_fig1; gca];
hold on;

% 1. Plot shaded uncertainty for the overall mean VSR
ci_lower_vsr = vsr_lin - vsr_lin_se;
ci_upper_vsr = vsr_lin + vsr_lin_se;
x_fill_vsr = [t_plot(:); flipud(t_plot(:))];
y_fill_vsr = [t_plot(:)*0 + ci_lower_vsr; t_plot(:)*0 + ci_upper_vsr];
fill(x_fill_vsr, y_fill_vsr, 0.5*[1 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2, ...
    'DisplayName', 'Overall Mean \pm SE');

% 2. Plot shaded uncertainty for the Linear Trend VSR
ci_lower_trend = vsr_trend - vsr_trend_se;
ci_upper_trend = vsr_trend + vsr_trend_se;
x_fill_trend = [t_plot(:); flipud(t_plot(:))];
y_fill_trend = [ci_lower_trend(:); flipud(ci_upper_trend(:))];
fill(x_fill_trend, y_fill_trend, [0.4 0.8 0.4], 'EdgeColor', 'none', 'FaceAlpha', 0.3, ...
    'DisplayName', 'Linear Fit \pm SE');

% 3. Plot overall mean VSR line
plot(t_plot, t_plot*0 + vsr_lin, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'Overall Mean VSR');

% 4. Plot the Linear Trend VSR
plot(t_plot, vsr_trend, 'g--', 'LineWidth', 2, ...
    'DisplayName', 'Linear Fit VSR');

% 5. Plot the end-point extrapolated VSR
plot(t_plot, vsr_interp, 'r-', 'LineWidth', 1.5, ...
    'DisplayName', 'Extrapolated VSR (End-point)');

% 6. Plot sliding window VSR estimates
errorbar(vvel.t, vvel.vsr_lin, vvel.vsr_lin_se, 'b.-', 'MarkerSize', 12, 'LineWidth', 1, ...
    'DisplayName', 'Linear VSR (Sliding)');

datetick('x', 'mm/yy', 'keeplimits');
ylabel('VSR (a^{-1})');
title('Linear Vertical Strain Rate over Time');
legend('Location', 'best');
grid on;

% Link x-axes for Figure 1
linkaxes(h_fig1, 'x');


%% 4. Basal Melt Rate Figure (4 Panels)
figure('Name', 'Basal Melt Rate Summary', 'Position', [150, 150, 800, 1000]);

% --- Panel 1: BMR with Constant VSR ---
% Calculation: BMR = -1 * (Bed Rate - Constant VSR)
bmr_const = -1 * (y_bed_merged_yr - vsr_lin);
bmr_const_err = sqrt(y_bed_ci_yr.^2 + vsr_lin_se^2);

subplot(4,1,1); h_fig2 = gca;
hold on;
errorbar(t_plot, bmr_const, bmr_const_err, 'Color', [0.8 0.8 0.8], 'LineStyle', 'none', ...
    'CapSize', 0, 'DisplayName', 'Geometric Uncertainty');
plot(t_plot, bmr_const, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, ...
    'DisplayName', 'BMR (Constant VSR)');
datetick('x', 'mm/yy', 'keeplimits');
ylabel('BMR (m/a)');
title('Basal Melt Rate (High-Res Bed, Constant Overall VSR)');
legend('Location', 'best');
grid on;

% --- Panel 2: BMR with Linear Fit VSR ---
% Calculation: BMR = -1 * (Bed Rate - Trend VSR)
bmr_trend = -1 * (y_bed_merged_yr - vsr_trend);
bmr_trend_err = sqrt(y_bed_ci_yr.^2 + vsr_trend_se.^2);

subplot(4,1,2); h_fig2 = [h_fig2; gca];
hold on;
errorbar(t_plot, bmr_trend, bmr_trend_err, 'Color', [0.8 0.8 0.8], 'LineStyle', 'none', ...
    'CapSize', 0, 'DisplayName', 'Geometric Uncertainty');
plot(t_plot, bmr_trend, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, ...
    'DisplayName', 'BMR (Linear Fit VSR)');
datetick('x', 'mm/yy', 'keeplimits');
ylabel('BMR (m/a)');
title('Basal Melt Rate (High-Res Bed, Linear Fit VSR)');
legend('Location', 'best');
grid on;

% --- Panel 3: BMR with End-Point Extrapolated Time-Variable VSR ---
% Calculation: BMR = -1 * (Bed Rate - Interpolated VSR)
bmr_interp = -1 * (y_bed_merged_yr - vsr_interp);
bmr_interp_err = sqrt(y_bed_ci_yr.^2 + vsr_se_interp.^2);

subplot(4,1,3); h_fig2 = [h_fig2; gca];
hold on;
errorbar(t_plot, bmr_interp, bmr_interp_err, 'Color', [0.8 0.8 0.8], 'LineStyle', 'none', ...
    'CapSize', 0, 'DisplayName', 'Geometric Uncertainty');
plot(t_plot, bmr_interp, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, ...
    'DisplayName', 'BMR (Interpolated/End-point VSR)');
datetick('x', 'mm/yy', 'keeplimits');
ylabel('BMR (m/a)');
title('Basal Melt Rate (High-Res Bed, End-Point Extrapolated VSR)');
legend('Location', 'best');
grid on;

% --- Panel 4: BMR strictly using Sliding Window Data ---
% Calculation: BMR = -1 * (Sliding Bed Rate - Sliding VSR)
bmr_slide = -1 * (vvel.v_bed - vvel.vsr_lin);
bmr_slide_err = sqrt(vvel.v_bed_se.^2 + vvel.vsr_lin_se.^2);

subplot(4,1,4); h_fig2 = [h_fig2; gca];
hold on;
errorbar(vvel.t, bmr_slide, bmr_slide_err, 'm.-', 'MarkerSize', 12, 'LineWidth', 1.5, ...
    'DisplayName', 'BMR (Sliding Window)');
datetick('x', 'mm/yy', 'keeplimits');
ylabel('BMR (m/a)');
title('Basal Melt Rate (Sliding Window Estimates)');
legend('Location', 'best');
grid on;

% Link x-axes for Figure 2
linkaxes(h_fig2, 'xy');