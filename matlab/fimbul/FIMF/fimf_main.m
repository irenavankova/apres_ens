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

% Sliding Window parameters (in days)
dt_window = 30; 
dt_over = 5;    

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
plot_vsr_results = true; 

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
    'plot_results', plot_vsr_results);

% Display overall mean results
fprintf('--- VSR Estimation Results at Bed Depth (z = %.2f m) ---\n', z_bed_comb);
fprintf('Bed Rate of Change: %.4f +/- %.4f m/a\n', v_bed, v_bed_se);
fprintf('Linear VSR thinning rate:         %.6f +/- %.6f m/a\n', vsr_lin, vsr_lin_se);
fprintf('Quadratic VSR thinning rate:      %.6f +/- %.6f m/a\n', vsr_quad, vsr_quad_se);

%% 2. Run VSR estimation over SLIDING WINDOW
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
if plot_vsr_results
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
        'plot_results', plot_vsr_results, ...
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
    
    if plot_vsr_results
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

% Setup plotting limits
t_plot = t_bed(1:length(y_bed_merged_yr));

% Subplot 1: Bed Rate of Change vs y_bed_merged
subplot(2,1,1);
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
subplot(2,1,2);
hold on;

% 1. Plot shaded uncertainty for the overall mean VSR
ci_lower_vsr = vsr_lin - vsr_lin_se;
ci_upper_vsr = vsr_lin + vsr_lin_se;
x_fill_vsr = [t_plot(:); flipud(t_plot(:))];
y_fill_vsr = [t_plot(:)*0 + ci_lower_vsr; t_plot(:)*0 + ci_upper_vsr];

fill(x_fill_vsr, y_fill_vsr, 0.5*[1 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2, ...
    'DisplayName', 'Overall Mean \pm SE');

% 2. Plot overall mean VSR line
plot(t_plot, t_plot*0 + vsr_lin, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'Overall Mean VSR');

% 3. Plot sliding window VSR estimates
errorbar(vvel.t, vvel.vsr_lin, vvel.vsr_lin_se, 'b.-', 'MarkerSize', 12, 'LineWidth', 1, ...
    'DisplayName', 'Linear VSR (Sliding)');

datetick('x', 'mm/yy', 'keeplimits');
ylabel('VSR (a^{-1})');
title('Linear Vertical Strain Rate over Time');
legend('Location', 'best');
grid on;
