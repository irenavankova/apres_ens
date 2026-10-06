% FIMF_LOAD_AND_PLOT
% Reads the saved variables from the NetCDF file and generates final figures.

% Set the filename corresponding to the site run
sitename = 'FIMFY25_B';
nc_filename = sprintf('%s_BMR_results.nc', sitename);

% Check if file exists
if ~isfile(nc_filename)
    error('File %s not found. Please run fimf_main first.', nc_filename);
end

% Load Variables
t_bed           = ncread(nc_filename, 't_bed');
y_bed_merged_yr = ncread(nc_filename, 'y_bed_merged_yr');
y_bed_ci_yr     = ncread(nc_filename, 'y_bed_ci_yr');
vsr_lin         = ncread(nc_filename, 'vsr_lin');
vsr_lin_se      = ncread(nc_filename, 'vsr_lin_se');

vvel_t          = ncread(nc_filename, 'vvel_t');
vvel_vsr_lin    = ncread(nc_filename, 'vvel_vsr_lin');
vvel_vsr_lin_se = ncread(nc_filename, 'vvel_vsr_lin_se');
vvel_v_bed      = ncread(nc_filename, 'vvel_v_bed');
vvel_v_bed_se   = ncread(nc_filename, 'vvel_v_bed_se');

% Call Plotting Function
fimf_plot_final_results(t_bed, y_bed_merged_yr, y_bed_ci_yr, vsr_lin, vsr_lin_se, ...
    vvel_t, vvel_vsr_lin, vvel_vsr_lin_se, vvel_v_bed, vvel_v_bed_se);

disp('Figures successfully generated.');