%% Main script: Multi-frequency VSR mean estimation

sitename = 'FIMFY25_B';
f0 = 200;
fend = 360;
db = 10;
bw = 40;
t1 = []; 
t2 = []; 
drnf_vsr = [0 40 ;124 4000]; % Depth ranges to exclude
drnf_artefacts = [];

% USER OPTIONS
opt_combine_method = 'median'; % 'mean' or 'median' for combining BW estimates
opt_remove_outliers_vel = true; % Remove outliers from velocity estimates across BWs before combining
opt_vel_method = 'line_fit'; % 'deriv_outlier', 'deriv_simple', or 'line_fit'
opt_deriv_combine = 'median'; % 'mean' or 'median' if opt_vel_method is 'deriv_simple'
opt_bed_source = 'tp'; % 'tp' (loads bed_tpq) or 'xcor' (loads bed_xcor)

% PLOTTING PARAMETERS
ylim_val = [0 170]; % Adjustable vertical y-axis limits (Range / Depth)
ytick_dist = 20;    % Adjustable distance between y-axis ticks

f1 = f0:db:fend;
N = length(f1);

% Matrices to accumulate frequency instance data
all_int_vel = [];
all_bed_vel = [];
all_bed_z = [];
all_ampCor_mean = [];
all_ampCor_std = [];

% Placeholder for external NaN filter array (e.g., from fimf_analysis)
iNaN_array = []; 

for k = 1:N
    sname = [sitename num2str(f1(k)) num2str(f1(k)+bw)];
    disp(['Processing: ' sname]);
    
    % Strict file loading for internals
    mat1 = [sname '_ts_fine_TT.mat'];
    if ~exist(mat1, 'file')
        error('File %s not found. Cannot proceed.', mat1);
    end
    load(mat1, 'ct');
    
    % Strict file loading for basal reflector based on user choice
    if strcmp(opt_bed_source, 'tp')
        mat2 = [sname '_bed_tpq.mat'];
        if ~exist(mat2, 'file')
            error('File %s not found. Cannot proceed.', mat2);
        end
        bed_data = load(mat2, 'tp');
        bed_struct = bed_data.tp;
    else
        mat2 = [sname '_bed_xcor.mat'];
        if ~exist(mat2, 'file')
            error('File %s not found. Cannot proceed.', mat2);
        end
        bed_data = load(mat2, 'bed_xcor');
        bed_struct = bed_data.bed_xcor;
    end
    
    % Setup Time Bounds
    time = ct.dts_xcor.time;
    if isempty(t1), t1_use = time(1); else, t1_use = t1; end
    if isempty(t2), t2_use = time(end); else, t2_use = t2; end
    tind = find(time >= t1_use & time <= t2_use);
    
    % Apply provided filtering (NaN exclusion)
    if ~isempty(iNaN_array)
        invalid_idx = intersect(tind, iNaN_array); 
        ct.dts_xcor.dh(invalid_idx, :) = NaN;
        if strcmp(opt_bed_source, 'tp')
            bed_struct.thickness(invalid_idx) = NaN;
        else
            bed_struct.dh(invalid_idx) = NaN;
        end
    end
    
    % Store depth arrays 
    z_all = ct.dts_xcor.dhRange;
    if strcmp(opt_bed_source, 'tp')
        all_bed_z = cat(1, all_bed_z, bed_struct.x_input);
    else
        all_bed_z = cat(1, all_bed_z, bed_struct.dhRange);
    end
    
    % Calculate Mean Velocities for this frequency instance
    [v_int, v_bed] = fimf_calc_dhdt(ct.dts_xcor, bed_struct, tind, opt_vel_method, opt_deriv_combine, opt_bed_source);
    
    % Accumulate results
    all_int_vel = cat(1, all_int_vel, v_int);
    all_bed_vel = cat(1, all_bed_vel, v_bed);
    
    % Ensure ampCor indices do not exceed its N-1 size
    tind_amp = tind(tind <= size(ct.dts_xcor.ampCor, 1));
    all_ampCor_mean = cat(1, all_ampCor_mean, mean(ct.dts_xcor.ampCor(tind_amp, :), 1, 'omitnan'));
    all_ampCor_std = cat(1, all_ampCor_std, std(ct.dts_xcor.ampCor(tind_amp, :), 1, 'omitnan'));
end

%% Statistics across central frequencies
% Remove outliers before taking the representative velocity and ranges
if opt_remove_outliers_vel
    ind_out_int = isoutlier(all_int_vel, 'median', 1);
    all_int_vel(ind_out_int) = NaN;
    
    ind_out_bed = isoutlier(all_bed_vel, 'median', 1);
    all_bed_vel(ind_out_bed) = NaN;
    
    ind_out_bed_z = isoutlier(all_bed_z, 'median', 1);
    all_bed_z(ind_out_bed_z) = NaN;
end

% Apply user specified combine method for VSR estimation and bed positioning
if strcmp(opt_combine_method, 'median')
    v_all = median(all_int_vel, 1, 'omitnan');
    v_bed_comb = median(all_bed_vel, 1, 'omitnan');
    z_bed_comb = median(all_bed_z, 1, 'omitnan');
else
    v_all = mean(all_int_vel, 1, 'omitnan');
    v_bed_comb = mean(all_bed_vel, 1, 'omitnan');
    z_bed_comb = mean(all_bed_z, 1, 'omitnan');
end

% Use standard deviation divided by sqrt(N) for ensemble error weights
v_all_se = std(all_int_vel, 0, 1, 'omitnan') ./ sqrt(N);

mean_ampCor = mean(all_ampCor_mean, 1, 'omitnan');
std_ampCor = mean(all_ampCor_std, 1, 'omitnan');

%% VSR Calculation and Fitting
[ind_fit, q_best, q_se, p_best, p_best_se] = fimf_calc_vsr(z_all', v_all', v_all_se', drnf_vsr);

%% Final Plot Generation
fimf_plot_vsr(sitename, z_all, all_bed_z, v_all, all_int_vel, v_bed_comb, z_bed_comb, all_bed_vel, v_all_se, ind_fit, q_best, q_se, p_best, p_best_se, all_ampCor_mean, mean_ampCor, std_ampCor, drnf_artefacts, ylim_val, ytick_dist);