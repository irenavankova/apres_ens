%Taken from TI05_S_vsr_estimation
%% Load site and calculate mean velocites

sname = 'FIMFY25';
opt_print = 0;

%spec = fs_site_specs(sname);

spec.t1 = [];
spec.t2 = [];
spec.drnf_vsr = [120 4000];
spec.depthEdges = [128 136];
spec.yrs = {'Y25'};
spec.vsr_removal = 'lin_all';
spec.mrLP_joint = 1;

drnf_vsr = spec.drnf_vsr;
drnf_artefacts = [];

try 
    load([sname '_ball.mat'])
    st = site;
end
load([sname '.mat'])

load([sname '_ts_fine_TT.mat'])
%ct_plot_xcor_vs_uwrp_lines(ct)

time = ct.dts_xcor.time;

t1 = spec.t1;
t2 = spec.t2;
if isempty(t2)
    t2 = time(end);
end
if isempty(t1)
    t1 = time(1);
end
tind = find(time >= t1 & time <= t2);
tl = [min(tind) max(tind)];

[hax,v_all,z_all,v_bed,z_bed,p_best,q_best,ind_fit,p_best_se,v_all_se] = A1_plots_vsr([],ct,[],[],[],drnf_vsr,0,[],0);

%% Final plot
[hax,plotwidth] = sub_sub(4,1,3,10,0.25,0);
cmap_jet = colormap(jet(64));

set(gcf,'CurrentAxes',hax(1))
%load(['TI05Y17_ball_c30.mat']);
amp = 20*log10(abs(site.spec_cor(:,:)'));
plot(mean(amp,2),site.cfg.Rcoarse,'k'); hold on
try
    amp = 20*log10(abs(st.spec_cor(:,:)'));
    plot(mean(amp,2),st.cfg.Rcoarse,'r'); hold on
    %legend('spc mean','burst mean')
end
xlim([-150 -40])
xlabel('Amplitude (dB)')

set(gcf,'CurrentAxes',hax(2))
mm = ct.dts_xcor.dh;
mstd = std(mm,[],1);
plot(mstd,ct.dts_xcor.dhRange, '-','color','k'); hold on
xlabel('std dev displ (m)')
xlim([0 2.8])

set(gcf,'CurrentAxes',hax(3))
plot(mean(ct.dts_xcor.ampCor),ct.dts_xcor.dhRange,'k','LineWidth',1); hold on
patch([mean(ct.dts_xcor.ampCor)-std(ct.dts_xcor.ampCor) fliplr(mean(ct.dts_xcor.ampCor)+std(ct.dts_xcor.ampCor))],[ct.dts_xcor.dhRange fliplr(ct.dts_xcor.dhRange)],'k','Facealpha',0.2,'Edgecolor','none'); hold on
xlim([0.75 1.01])
xlabel('Correlation coef')
for j = 4%:5
    set(gcf,'CurrentAxes',hax(j))
    
    w = 1./v_all_se;
    %w = 1./mean(ct.dts_xcor.ampCor);
    w = w/max(w);
    [ind,q_best,q_se] = ct_tide_quadratic_fit(z_all,v_all,w,drnf_vsr,'lsq_weighted');
    [~,p_best,~,p_best_se] = ct_tide_line_fit(z_all,v_all,w,drnf_vsr,'lsq_weighted');

    plot(v_all,z_all,'.','color',cmap_jet(44,:)); hold on
    plot(v_all(ind),z_all(ind),'.','color',cmap_jet(6,:)); hold on

    q = q_best+ q_se;
    yQL = ct_tide_compute_quadratic(z_all,q); % Best quad fit
    q = q_best-q_se;
    yQU = ct_tide_compute_quadratic(z_all,q); % Best quad fit
    yQ = ct_tide_compute_quadratic(z_all,q_best); % Best quad fit
    patch([yQL fliplr(yQU)],[z_all fliplr(z_all)],'k','Facealpha',0.2,'Edgecolor','none'); hold on
    plot(yQ,z_all,'k','LineWidth',1); hold on

    yL = p_best(1)*z_all + p_best(2); % Best line fit
    pU = p_best - p_best_se;
    yLU = pU(1)*z_all + pU(2); % Best line fit
    pL = p_best + p_best_se;
    yLL = pL(1)*z_all + pL(2); % Best line fit
    patch([yLL fliplr(yLU)],[z_all fliplr(z_all)],'r','Facealpha',0.2,'Edgecolor','none'); hold on
    plot(yL,z_all,'r','LineWidth',1); hold on

    xlim([-2.5 6])    
    xlabel('vel (m a^{-1})')
    % Calculate Ftest values
    y_data = v_all(ind);
    y_lin = yL(ind);
    y_quad = yQ(ind);
    %y_data = y_lin;
    [dof,faj,ftab] = ct_calc_Ftest_QL_nested(y_data,y_lin,y_quad);

    %w = w*0+1;
    [dof2,faj2,ftab2] = tg_rev_Ftest_weighted(y_data,y_lin,y_quad,w(ind));

    [[dof faj ftab]' [dof2 faj2 ftab2]']
    text(0.6,0.96,['F = ' num2str(round(faj2,1))],'Units','normalized','fontsize',8,'FontWeight','normal')

end
for j = 1:length(hax)
    set(gcf,'CurrentAxes',hax(j))
    set(gca,'Fontsize',8,'ydir','reverse')
    ct_tide_plot_shade_drnf(xlim,drnf_artefacts,'b')
end
abcd = 'abcde';
for j = 1:length(hax)
    set(gcf,'CurrentAxes',hax(j))
    set(gca,'Fontsize',8,'Linewidth',0.5,'ytick',0:100:4000,'ydir','reverse')
    if j > 1
        set(gca,'YTickLabel',[])
        ylabel('')
    else
        ylabel('Range (m)')
    end
    text(0.05,0.96,[abcd(j)],'Units','normalized','fontsize',8,'FontWeight','normal')
    axis tight
    %ylim([0 1096])
end
linkaxes(hax,'y')

xlim([-2 1.5])
if opt_print
    print(gcf,'-r300', '-dpng',[ sname '_vsr_mean_estimation_Aall.png']);
end

