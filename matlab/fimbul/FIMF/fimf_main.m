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
opt_bed_source = 'xcor'; % 'tp' (loads bed_tpq) or 'xcor' (loads bed_xcor)

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
drnf_vsr = [0 20; 124 4000]; %[0 40; 124 4000]
drnf_artefacts = [];
max_nan_threshold = 0.5; % More than 20% bad -> do NaN

% Sliding Window parameters (in days)
dt_window = 30; 
dt_over = 1;

% Output and saving options
opt_save_to_nc = false;

% Load all raw data once
[bed_data, int_data, t_bed, t_int, f1] = fimf_load_matfiles(...
    sitename, f0, fend, db, bw, i1t, dhRange_low, dhRange_high, opt_bed_source);

% Run filter jumps in 'final' mode to get iNaN_array and uncertainties
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
        'hax', hax_anim); 
        
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

%% 3. Format Data and Save to NetCDF
% Convert arrays to consistent units and shapes for saving
y_bed_merged_yr = y_bed_merged(:) * d2y;
y_bed_ci_yr = y_bed_ci(:) * d2y;
t_bed_out = t_bed(1:end-1); % Match length to difference vector
t_bed_out = t_bed_out(:);

% Sliding window vectors
vvel_t = vvel.t(:);
vvel_vsr_lin = vvel.vsr_lin(:);
vvel_vsr_lin_se = vvel.vsr_lin_se(:);
vvel_v_bed = vvel.v_bed(:);
vvel_v_bed_se = vvel.v_bed_se(:);

if opt_save_to_nc == true
    disp('--- Saving Results to NetCDF ---');
    
    % File setup
    nc_filename = sprintf('%s_BMR_results.nc', sitename);
    if isfile(nc_filename)
        delete(nc_filename);
    end
    
    % Define dimensions
    dim_bed = length(t_bed_out);
    dim_slide = length(vvel_t);
    
    % --- Create and write variables with descriptions ---
    
    % 1. t_bed
    nccreate(nc_filename, 't_bed', 'Dimensions', {'time_bed', dim_bed});
    ncwrite(nc_filename, 't_bed', t_bed_out);
    ncwriteatt(nc_filename, 't_bed', 'description', 'Time vector corresponding to high-resolution bed measurements');
    ncwriteatt(nc_filename, 't_bed', 'units', 'datenum / days');
    
    % 2. y_bed_merged_yr
    nccreate(nc_filename, 'y_bed_merged_yr', 'Dimensions', {'time_bed', dim_bed});
    ncwrite(nc_filename, 'y_bed_merged_yr', y_bed_merged_yr);
    ncwriteatt(nc_filename, 'y_bed_merged_yr', 'description', 'Bed range (total thickness) rate of change timeseries');
    ncwriteatt(nc_filename, 'y_bed_merged_yr', 'units', 'm/a');
    
    % 3. y_bed_ci_yr
    nccreate(nc_filename, 'y_bed_ci_yr', 'Dimensions', {'time_bed', dim_bed});
    ncwrite(nc_filename, 'y_bed_ci_yr', y_bed_ci_yr);
    ncwriteatt(nc_filename, 'y_bed_ci_yr', 'description', '95% Confidence interval for bed range (total thickness) rate of change');
    ncwriteatt(nc_filename, 'y_bed_ci_yr', 'units', 'm/a');
    
    % 4. vsr_lin
    nccreate(nc_filename, 'vsr_lin', 'Dimensions', {'scalar', 1});
    ncwrite(nc_filename, 'vsr_lin', vsr_lin);
    ncwriteatt(nc_filename, 'vsr_lin', 'description', 'Time average strain thinning rate (assuming constant vertical strain rate profile) estimated over the whole time series');
    ncwriteatt(nc_filename, 'vsr_lin', 'units', 'a^-1');
    
    % 5. vsr_lin_se
    nccreate(nc_filename, 'vsr_lin_se', 'Dimensions', {'scalar', 1});
    ncwrite(nc_filename, 'vsr_lin_se', vsr_lin_se);
    ncwriteatt(nc_filename, 'vsr_lin_se', 'description', 'Standard error of time average strain thinning rate (assuming constant vertical strain rate profile) over the whole time series');
    ncwriteatt(nc_filename, 'vsr_lin_se', 'units', 'a^-1');
    
    % 6. vvel_t
    nccreate(nc_filename, 'vvel_t', 'Dimensions', {'time_slide', dim_slide});
    ncwrite(nc_filename, 'vvel_t', vvel_t);
    ncwriteatt(nc_filename, 'vvel_t', 'description', 'Mid-window timestamps for sliding window strain thinning rate timeseries estimates');
    ncwriteatt(nc_filename, 'vvel_t', 'units', 'datenum / days');
    
    % 7. vvel_vsr_lin
    nccreate(nc_filename, 'vvel_vsr_lin', 'Dimensions', {'time_slide', dim_slide});
    ncwrite(nc_filename, 'vvel_vsr_lin', vvel_vsr_lin);
    ncwriteatt(nc_filename, 'vvel_vsr_lin', 'description', 'Strain thinning rate (assuming constant vertical strain rate profile) calculated over sliding window');
    ncwriteatt(nc_filename, 'vvel_vsr_lin', 'units', 'a^-1');
    
    % 8. vvel_vsr_lin_se
    nccreate(nc_filename, 'vvel_vsr_lin_se', 'Dimensions', {'time_slide', dim_slide});
    ncwrite(nc_filename, 'vvel_vsr_lin_se', vvel_vsr_lin_se);
    ncwriteatt(nc_filename, 'vvel_vsr_lin_se', 'description', 'Standard error of sliding window strain thinning rate estimate');
    ncwriteatt(nc_filename, 'vvel_vsr_lin_se', 'units', 'a^-1');
    
    % 9. vvel_v_bed
    nccreate(nc_filename, 'vvel_v_bed', 'Dimensions', {'time_slide', dim_slide});
    ncwrite(nc_filename, 'vvel_v_bed', vvel_v_bed);
    ncwriteatt(nc_filename, 'vvel_v_bed', 'description', 'Bed range rate of change estimated over sliding window');
    ncwriteatt(nc_filename, 'vvel_v_bed', 'units', 'm/a');
    
    % 10. vvel_v_bed_se
    nccreate(nc_filename, 'vvel_v_bed_se', 'Dimensions', {'time_slide', dim_slide});
    ncwrite(nc_filename, 'vvel_v_bed_se', vvel_v_bed_se);
    ncwriteatt(nc_filename, 'vvel_v_bed_se', 'description', 'Standard error of sliding window bed range rate of change');
    ncwriteatt(nc_filename, 'vvel_v_bed_se', 'units', 'm/a');
    
    % Optional: Global attributes for overall metadata
    ncwriteatt(nc_filename, '/', 'site_name', sitename);
    ncwriteatt(nc_filename, '/', 'creation_date', datestr(now));
    
    disp(['Results successfully saved to ', nc_filename]);

end


%% 4. Generate Final Plots
fimf_plot_final_results(t_bed_out, y_bed_merged_yr, y_bed_ci_yr, vsr_lin, vsr_lin_se, ...
                        vvel_t, vvel_vsr_lin, vvel_vsr_lin_se, vvel_v_bed, vvel_v_bed_se);