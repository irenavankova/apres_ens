%% Main script: dh/dt and its variance across bandwidth
% Parameters
bw = 40;
db = 10;
f0 = 200;
fend = 360;
Nbad = 6;
Nbaddh = 3;
i1t = 3;
sitename = 'FIMFY25_B';
dhRange_low = 20;
dhRange_high = 80;
plot_results = false;
ymax = 0.06;
Lwide = 1;

% Define parameters
num_bins = 20;
quant_grid = linspace(0.05, 0.95, 20);

% Load data for bed and internal cases
[t_bed, dhdtmat_bed, ~] = fimf_load_and_create_timeseries(sitename, f0, fend, db, bw, i1t, true, dhRange_low, dhRange_high);
[t_int, dhdtmat_int, ~] = fimf_load_and_create_timeseries(sitename, f0, fend, db, bw, i1t, false, dhRange_low, dhRange_high);

% Compute raw, std, and confidence intervals in the main script
y_bed_raw = mean(dhdtmat_bed, 2, 'omitnan');
y_bed_std = std(dhdtmat_bed, [], 2, 'omitnan');
y_int_raw = mean(dhdtmat_int, 2, 'omitnan');
y_int_std = std(dhdtmat_int, [], 2, 'omitnan');

[~, ~, ci_bed] = ttest(dhdtmat_bed', 0, 'Alpha', 0.05); 
y_bed_ci_lower = ci_bed(1, :)'; y_bed_ci_upper = ci_bed(2, :)';
y_bed_ci = (y_bed_ci_upper - y_bed_ci_lower) / 2;

[~, ~, ci_int] = ttest(dhdtmat_int', 0, 'Alpha', 0.05); 
y_int_ci_lower = ci_int(1, :)'; y_int_ci_upper = ci_int(2, :)';
y_int_ci = (y_int_ci_upper - y_int_ci_lower) / 2;

%% Bandwidth filter
bw_med_filt_param = 5;
[iNaN_bw, iNaN_bed_bw, iNaN_int_bw] = fimf_filter_bw(dhdtmat_bed, dhdtmat_int, t_bed, t_int, Nbad, bw_med_filt_param, y_bed_raw, y_bed_std, y_int_raw, y_int_std, ymax, plot_results);

y_bed_bw = y_bed_raw; y_bed_bw(iNaN_bw) = NaN;
y_bed_bw_bed = y_bed_raw; y_bed_bw_bed(iNaN_bed_bw) = NaN;
y_bed_bw_int = y_bed_raw; y_bed_bw_int(iNaN_int_bw) = NaN;

y_int_bw = y_int_raw; y_int_bw(iNaN_bw) = NaN;
y_int_bw_bed = y_int_raw; y_int_bw_bed(iNaN_bed_bw) = NaN;
y_int_bw_int = y_int_raw; y_int_bw_int(iNaN_int_bw) = NaN;

%% Moving median filter
window_size = 31 * 13;
[iNaN_movmedian, iNaN_bed_movmedian, iNaN_int_movmedian] = fimf_filter_movmedian(dhdtmat_bed, dhdtmat_int, t_bed, t_int, window_size, y_bed_raw, y_bed_std, y_int_raw, y_int_std, ymax, plot_results);

y_bed_movmedian = y_bed_raw; y_bed_movmedian(iNaN_movmedian) = NaN;
y_bed_movmedian_bed = y_bed_raw; y_bed_movmedian_bed(iNaN_bed_movmedian) = NaN;
y_bed_movmedian_int = y_bed_raw; y_bed_movmedian_int(iNaN_int_movmedian) = NaN;

y_int_movmedian = y_int_raw; y_int_movmedian(iNaN_movmedian) = NaN;
y_int_movmedian_bed = y_int_raw; y_int_movmedian_bed(iNaN_bed_movmedian) = NaN;
y_int_movmedian_int = y_int_raw; y_int_movmedian_int(iNaN_int_movmedian) = NaN;

%% Median filter
[iNaN_median, iNaN_bed_median, iNaN_int_median] = fimf_filter_median(dhdtmat_bed, dhdtmat_int, t_bed, t_int, y_bed_raw, y_bed_std, y_int_raw, y_int_std, ymax, plot_results);

y_bed_median = y_bed_raw; y_bed_median(iNaN_median) = NaN;
y_bed_median_bed = y_bed_raw; y_bed_median_bed(iNaN_bed_median) = NaN;
y_bed_median_int = y_bed_raw; y_bed_median_int(iNaN_int_median) = NaN;

y_int_median = y_int_raw; y_int_median(iNaN_median) = NaN;
y_int_median_bed = y_int_raw; y_int_median_bed(iNaN_bed_median) = NaN;
y_int_median_int = y_int_raw; y_int_median_int(iNaN_int_median) = NaN;

%% KS filter
ref_len_KS = floor(length(y_bed_raw)/2);
window_len_KS = floor(ref_len_KS/10);
[iNaN_ks, iNaN_bed_ks, iNaN_int_ks] = fimf_filter_ks(...
    dhdtmat_bed, dhdtmat_int, t_bed, t_int, y_bed_raw, y_bed_std, y_int_raw, y_int_std, ymax, window_len_KS, ref_len_KS, plot_results);

y_bed_ks = y_bed_raw; y_bed_ks(iNaN_ks) = NaN;
y_bed_ks_bed = y_bed_raw; y_bed_ks_bed(iNaN_bed_ks) = NaN;
y_bed_ks_int = y_bed_raw; y_bed_ks_int(iNaN_int_ks) = NaN;

y_int_ks = y_int_raw; y_int_ks(iNaN_ks) = NaN;
y_int_ks_bed = y_int_raw; y_int_ks_bed(iNaN_bed_ks) = NaN;
y_int_ks_int = y_int_raw; y_int_ks_int(iNaN_int_ks) = NaN;

%% KS filter (Fixed Reference)
ref_start = 500; 
ref_len_KS_fix = floor(length(y_bed_raw)/20);
window_len_KS_fix = floor(ref_len_KS_fix/5);
[iNaN_ks_fix, iNaN_bed_ks_fix, iNaN_int_ks_fix] = fimf_filter_ks_fix_ref(...
    dhdtmat_bed, dhdtmat_int, t_bed, t_int, y_bed_raw, y_bed_std, y_int_raw, y_int_std, ymax, window_len_KS_fix, ref_start, ref_len_KS_fix, plot_results);

y_bed_ks_fix = y_bed_raw; if ~isempty(iNaN_ks_fix), y_bed_ks_fix(iNaN_ks_fix) = NaN; end
y_bed_ks_fix_bed = y_bed_raw; y_bed_ks_fix_bed(iNaN_bed_ks_fix) = NaN;
y_bed_ks_fix_int = y_bed_raw; y_bed_ks_fix_int(iNaN_int_ks_fix) = NaN;

y_int_ks_fix = y_int_raw; if ~isempty(iNaN_ks_fix), y_int_ks_fix(iNaN_ks_fix) = NaN; end
y_int_ks_fix_bed = y_int_raw; y_int_ks_fix_bed(iNaN_bed_ks_fix) = NaN;
y_int_ks_fix_int = y_int_raw; y_int_ks_fix_int(iNaN_int_ks_fix) = NaN;

%% EMD filter (Sliding reference)
ref_len_EMD = floor(length(y_bed_raw)/50);
window_len_EMD = floor(ref_len_EMD/10);
[iNaN_emd, iNaN_bed_emd, iNaN_int_emd] = fimf_filter_emd(...
    dhdtmat_bed, dhdtmat_int, t_bed, t_int, y_bed_raw, y_bed_std, y_int_raw, y_int_std, ymax, window_len_EMD, ref_len_EMD, quant_grid, plot_results);

y_bed_emd = y_bed_raw; y_bed_emd(iNaN_emd) = NaN;
y_bed_emd_bed = y_bed_raw; y_bed_emd_bed(iNaN_bed_emd) = NaN;
y_bed_emd_int = y_bed_raw; y_bed_emd_int(iNaN_int_emd) = NaN;

y_int_emd = y_int_raw; y_int_emd(iNaN_emd) = NaN;
y_int_emd_bed = y_int_raw; y_int_emd_bed(iNaN_bed_emd) = NaN;
y_int_emd_int = y_int_raw; y_int_emd_int(iNaN_int_emd) = NaN;

%% EMD filter (Fixed Reference)
ref_start = 500;
ref_len_EMD_fix = floor(length(y_bed_raw)/50);
window_len_EMD_fix = floor(ref_len_EMD_fix/10);
[iNaN_emd_fix, iNaN_bed_emd_fix, iNaN_int_emd_fix] = fimf_filter_emd_fix_ref(...
    dhdtmat_bed, dhdtmat_int, t_bed, t_int, y_bed_raw, y_bed_std, y_int_raw, y_int_std, ymax, window_len_EMD_fix, ref_start, ref_len_EMD_fix, quant_grid, plot_results);

y_bed_emd_fix = y_bed_raw; if ~isempty(iNaN_emd_fix), y_bed_emd_fix(iNaN_emd_fix) = NaN; end
y_bed_emd_fix_bed = y_bed_raw; y_bed_emd_fix_bed(iNaN_bed_emd_fix) = NaN;
y_bed_emd_fix_int = y_bed_raw; y_bed_emd_fix_int(iNaN_int_emd_fix) = NaN;

y_int_emd_fix = y_int_raw; if ~isempty(iNaN_emd_fix), y_int_emd_fix(iNaN_emd_fix) = NaN; end
y_int_emd_fix_bed = y_int_raw; y_int_emd_fix_bed(iNaN_bed_emd_fix) = NaN;
y_int_emd_fix_int = y_int_raw; y_int_emd_fix_int(iNaN_int_emd_fix) = NaN;

%% JSD filter
ref_len_JSD = floor(length(y_bed_raw)/5);
window_len_JSD = floor(ref_len_JSD/30);
[iNaN_jsd, iNaN_bed_jsd, iNaN_int_jsd] = fimf_filter_jsd(...
    dhdtmat_bed, dhdtmat_int, t_bed, t_int, y_bed_raw, y_bed_std, y_int_raw, y_int_std, ymax, window_len_JSD, ref_len_JSD, num_bins, plot_results);

y_bed_jsd = y_bed_raw; y_bed_jsd(iNaN_jsd) = NaN;
y_bed_jsd_bed = y_bed_raw; y_bed_jsd_bed(iNaN_bed_jsd) = NaN;
y_bed_jsd_int = y_bed_raw; y_bed_jsd_int(iNaN_int_jsd) = NaN;

y_int_jsd = y_int_raw; y_int_jsd(iNaN_jsd) = NaN;
y_int_jsd_bed = y_int_raw; y_int_jsd_bed(iNaN_bed_jsd) = NaN;
y_int_jsd_int = y_int_raw; y_int_jsd_int(iNaN_int_jsd) = NaN;

%% JSD filter fixed reference window position
ref_len_JSD_fix = floor(length(y_bed_raw)/3);
window_len_JSD_fix = floor(ref_len_JSD_fix/80);
ref_start_JSD_fix = 500;
[iNaN_jsd_fix, iNaN_bed_jsd_fix, iNaN_int_jsd_fix] = fimf_filter_jsd_fix_ref(...
    dhdtmat_bed, dhdtmat_int, t_bed, t_int, y_bed_raw, y_bed_std, y_int_raw, y_int_std, ymax, window_len_JSD_fix, ref_start_JSD_fix, ref_len_JSD_fix, num_bins, plot_results);

y_bed_jsd_fix = y_bed_raw; y_bed_jsd_fix(iNaN_jsd_fix) = NaN;
y_bed_jsd_fix_bed = y_bed_raw; y_bed_jsd_fix_bed(iNaN_bed_jsd_fix) = NaN;
y_bed_jsd_fix_int = y_bed_raw; y_bed_jsd_fix_int(iNaN_int_jsd_fix) = NaN;

y_int_jsd_fix = y_int_raw; y_int_jsd_fix(iNaN_jsd_fix) = NaN;
y_int_jsd_fix_bed = y_int_raw; y_int_jsd_fix_bed(iNaN_bed_jsd_fix) = NaN;
y_int_jsd_fix_int = y_int_raw; y_int_jsd_fix_int(iNaN_int_jsd_fix) = NaN;

%% MAD filter
window_len_MAD = 91;
[iNaN_mad, iNaN_bed_mad, iNaN_int_mad] = fimf_filter_mad(...
    dhdtmat_bed, dhdtmat_int, t_bed, t_int, y_bed_raw, y_bed_std, y_int_raw, y_int_std, ymax, window_len_MAD, plot_results);

y_bed_mad = y_bed_raw; y_bed_mad(iNaN_mad) = NaN;
y_bed_mad_bed = y_bed_raw; y_bed_mad_bed(iNaN_bed_mad) = NaN;
y_bed_mad_int = y_bed_raw; y_bed_mad_int(iNaN_int_mad) = NaN;

y_int_mad = y_int_raw; y_int_mad(iNaN_mad) = NaN;
y_int_mad_bed = y_int_raw; y_int_mad_bed(iNaN_bed_mad) = NaN;
y_int_mad_int = y_int_raw; y_int_mad_int(iNaN_int_mad) = NaN;

%% --- COMPARISON FIGURES: All Filters for Bed and Internal Data ---
bed_filters = {
    y_bed_bw, 'BW';
    y_bed_movmedian, 'MovMed';
    y_bed_median, 'Median';
    y_bed_ks, 'KS';
    y_bed_emd, 'EMD';
    y_bed_jsd, 'JSD';
    y_bed_mad, 'MAD'
};
bed_filters_bed = {
    y_bed_bw_bed, 'BW';
    y_bed_movmedian_bed, 'MovMed';
    y_bed_median_bed, 'Median';
    y_bed_ks_bed, 'KS';
    y_bed_emd_bed, 'EMD';
    y_bed_jsd_bed, 'JSD';
    y_bed_mad_bed, 'MAD'
};
bed_filters_int = {
    y_bed_bw_int, 'BW';
    y_bed_movmedian_int, 'MovMed';
    y_bed_median_int, 'Median';
    y_bed_ks_int, 'KS';
    y_bed_emd_int, 'EMD';
    y_bed_jsd_int, 'JSD';
    y_bed_mad_int, 'MAD'
};

int_filters = {
    y_int_bw, 'BW';
    y_int_movmedian, 'MovMed';
    y_int_median, 'Median';
    y_int_ks, 'KS';
    y_int_emd, 'EMD';
    y_int_jsd, 'JSD';
    y_int_mad, 'MAD'
};
int_filters_bed = {
    y_int_bw_bed, 'BW';
    y_int_movmedian_bed, 'MovMed';
    y_int_median_bed, 'Median';
    y_int_ks_bed, 'KS';
    y_int_emd_bed, 'EMD';
    y_int_jsd_bed, 'JSD';
    y_int_mad_bed, 'MAD'
};
int_filters_int = {
    y_int_bw_int, 'BW';
    y_int_movmedian_int, 'MovMed';
    y_int_median_int, 'Median';
    y_int_ks_int, 'KS';
    y_int_emd_int, 'EMD';
    y_int_jsd_int, 'JSD';
    y_int_mad_int, 'MAD'
};

bed_filters_fix = {
    y_bed_bw, 'BW';
    y_bed_movmedian, 'MovMed';
    y_bed_median, 'Median';
    y_bed_emd_fix, 'EMDfix';
    y_bed_jsd_fix, 'JSDfix'
    };
bed_filters_fix_bed = {
    y_bed_bw_bed, 'BW';
    y_bed_movmedian_bed, 'MovMed';
    y_bed_median_bed, 'Median';
    y_bed_emd_fix_bed, 'EMDfix';
    y_bed_jsd_fix_bed, 'JSDfix'
    };
bed_filters_fix_int = {
    y_bed_bw_int, 'BW';
    y_bed_movmedian_int, 'MovMed';
    y_bed_median_int, 'Median';
    y_bed_emd_fix_int, 'EMDfix';
    y_bed_jsd_fix_int, 'JSDfix'
    };

int_filters_fix = {
    y_int_bw, 'BW';
    y_int_movmedian, 'MovMed';
    y_int_median, 'Median';
    y_int_emd_fix, 'EMDfix';
    y_int_jsd_fix, 'JSDfix'
    };
int_filters_fix_bed = {
    y_int_bw_bed, 'BW';
    y_int_movmedian_bed, 'MovMed';
    y_int_median_bed, 'Median';
    y_int_emd_fix_bed, 'EMDfix';
    y_int_jsd_fix_bed, 'JSDfix'
    };
int_filters_fix_int = {
    y_int_bw_int, 'BW';
    y_int_movmedian_int, 'MovMed';
    y_int_median_int, 'Median';
    y_int_emd_fix_int, 'EMDfix';
    y_int_jsd_fix_int, 'JSDfix'
    };

create_comparison_figure(bed_filters, t_bed, y_bed_raw, y_bed_std, y_bed_ci, ymax, 'Bed combo');
create_comparison_figure(bed_filters_bed, t_bed, y_bed_raw, y_bed_std, y_bed_ci, ymax, 'Bed bed');
create_comparison_figure(bed_filters_int, t_bed, y_bed_raw, y_bed_std, y_bed_ci, ymax, 'Bed int');
create_comparison_figure(int_filters, t_int, y_int_raw, y_int_std, y_int_ci, ymax, 'Int combo');
create_comparison_figure(int_filters_bed, t_int, y_int_raw, y_int_std, y_int_ci, ymax, 'Int combo');
create_comparison_figure(int_filters_int, t_int, y_int_raw, y_int_std, y_int_ci, ymax, 'Int combo');

create_comparison_figure(bed_filters_fix, t_bed, y_bed_raw, y_bed_std, y_bed_ci, ymax, 'BedGood combo');
create_comparison_figure(bed_filters_fix_bed, t_bed, y_bed_raw, y_bed_std, y_bed_ci, ymax, 'BedGood bed');
create_comparison_figure(bed_filters_fix_int, t_bed, y_bed_raw, y_bed_std, y_bed_ci, ymax, 'BedGood int');
create_comparison_figure(int_filters_fix, t_int, y_int_raw, y_int_ci, y_int_ci, ymax, 'IntGood combo');
create_comparison_figure(int_filters_fix_bed, t_int, y_int_raw, y_int_ci, y_int_ci, ymax, 'IntGood bed');
create_comparison_figure(int_filters_fix_int, t_int, y_int_raw, y_int_ci, y_int_ci, ymax, 'IntGood int');


%% Merge
bed_filters_final = {
    y_bed_bw, 'BW combo';
    y_bed_movmedian_int, 'MovMed int';
    y_bed_emd_fix_int, 'EMDfix int';
    y_bed_jsd_fix_int, 'JSDfix int'
    };

selected_series = cell2mat(bed_filters_final(:, 1)');
valid_count = sum(~isnan(selected_series), 2);
y_bed_merged = y_bed_raw;
y_bed_merged((size(selected_series, 2) - valid_count) >= 2) = NaN;
bed_filters_final = [bed_filters_final; {y_bed_merged, 'Merged'}];

int_filters_final = {
    y_int_bw, 'BW combo';
    y_int_movmedian_int, 'MovMed int';
    y_int_emd_fix_int, 'EMDfix int';
    y_int_jsd_fix_int, 'JSDfix int'
    };

selected_series_int = cell2mat(int_filters_final(:, 1)');
valid_count_int = sum(~isnan(selected_series_int), 2);
y_int_merged = y_int_raw;
y_int_merged((size(selected_series_int, 2) - valid_count_int) >= 2) = NaN;
int_filters_final = [int_filters_final; {y_int_merged, 'Merged'}];

create_comparison_figure(bed_filters_final, t_bed, y_bed_raw, y_bed_std, y_bed_ci, ymax, 'Bed final');
create_comparison_figure(int_filters_final, t_int, y_int_raw, y_int_std, y_int_ci, ymax, 'Int final');


%% Helper function for creating comparison figures with shaded CI
function create_comparison_figure(filters, t, y_raw, y_std, y_ci, ymax, title_str)
    num_filters = size(filters, 1);
    Lwide = 1;
    
    fig = figure;
    
    % Set figure units to centimeters
    set(fig, 'Units', 'centimeters');
    fig_pos = get(fig, 'Position');
    plotwidth = fig_pos(3);
    plotheight = fig_pos(4);
    
    % Define margins (in cm)
    leftmargin = 1.5;
    rightmargin = 1.0;
    bottommargin = 1.5;
    topmargin = 1.5;
    spacex = 0;
    spacey = 0.2; % 2 mm vertical spacing
    
    % Calculate normalized position vectors using subplot_pos
    pos = ap_subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, ...
                      bottommargin, topmargin, 1, num_filters, spacex, spacey);
                  
    h = zeros(num_filters, 1);
    t_plot = t(1:end-1);
    
    % Prepare vectors for polygon shading (lower bound forward, upper bound backward)
    ci_lower = y_raw - y_ci;
    ci_upper = y_raw + y_ci;
    x_fill = [t_plot(:); flipud(t_plot(:))];
    y_fill = [ci_lower(:); flipud(ci_upper(:))];
    
    for f = 1:num_filters
        j_idx = num_filters - f + 1;
        
        ax = axes('Position', pos{1, j_idx});
        h(f) = ax;
        
        y = filters{f, 1};
        label = filters{f, 2};
        
        % 1. Plot Light Gray Shaded Confidence Interval Region
        fill(x_fill, y_fill, 0.5*[1 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5); hold on;
        
        % 2. Plot Standard Deviation Bounds (Cyan Lines)
        %plot(t_plot, y_raw - y_std, 'color', 'cyan', 'LineWidth', Lwide);
        %plot(t_plot, y_raw + y_std, 'color', 'cyan', 'LineWidth', Lwide);
        
        % 3. Plot Filtered Data (Solid Black Line)
        plot(t_plot, y, 'k-', 'LineWidth', Lwide);
        
        grid on;
        ylabel(label, 'FontWeight', 'normal');
        datetick('x', 'mm/yy', 'keeplimits');
        ylim([-ymax ymax]);
        
        if f == 1
            title(title_str, 'fontweight', 'normal');
        end
        
        % Remove X-tick labels for all subplots except the bottom one
        if f < num_filters
            set(ax, 'XTickLabel', []);
        else
            xlabel('Time');
        end
    end
    
    linkaxes(h, 'xy');
end