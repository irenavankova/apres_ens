function [bed_data, int_data, t_bed, t_int, f1] = fimf_load_matfiles(sitename, f0, fend, db, bw, i1t, dhRange_low, dhRange_high, opt_bed_source)
% FIMF_LOAD_ALL_RAW_DATA Load all .mat files for bed and internal data.
%   [bed_data, int_data, t_bed, t_int, f1] = FIMF_LOAD_ALL_RAW_DATA(...)
%   - bed_data: Cell array of structs (tp or bed_xcor) for each frequency.
%   - int_data: Cell array of structs (ct) for each frequency.
%   - t_bed, t_int: Time vectors (column vectors).
%   - f1: Frequency vector.

f1 = f0:db:fend;
N = length(f1);

% Initialize outputs
bed_data = cell(1, N);
int_data = cell(1, N);
t_bed = [];
t_int = [];

for k = 1:N
    sname = [sitename num2str(f1(k)) num2str(f1(k)+bw)];

    % Load bed data
    if strcmp(opt_bed_source, 'tp')
        bed_file = [sname '_bed_tpq.mat'];
        bed_tmp = load(bed_file, 'tp');
        bed_data{k} = bed_tmp.tp;
        if isempty(t_bed)
            t_bed = bed_data{k}.tim(i1t:end);
            t_bed = t_bed(:);
        end
    else
        bed_file = [sname '_bed_xcor.mat'];
        bed_tmp = load(bed_file, 'bed_xcor');
        bed_file_tp = [sname '_bed_tpq.mat'];
        bed_tmp_tp = load(bed_file_tp, 'tp');
        bed_tmp.bed_xcor.x_input = bed_tmp_tp.tp.x_input;
        bed_data{k} = bed_tmp.bed_xcor;
        if isempty(t_bed)
            t_bed = bed_data{k}.time(i1t:end);
            t_bed = t_bed(:);
        end
    end

    % Load internal data
    int_file = [sname '_ts_fine_TT.mat'];
    int_tmp = load(int_file, 'ct');
    int_data{k} = int_tmp.ct;
    if isempty(t_int)
        t_int = int_data{k}.dts_xcor.time(i1t:end);
        t_int = t_int(:);
    end
end
end