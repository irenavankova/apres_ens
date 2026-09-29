function [hax,plotwidth,fontsize,f,plotheight] = ap_sub_sub(subplotsx,subplotsy,Lx,Ly,spacex,spacey,rightedge,topedge,leftedge,bottomedge)

fontsize=8; 

clear f

% Lx = 2.5;
% Ly = 10;
% 
% subplotsx=8;
% subplotsy=1;   

if nargin < 7
    leftedge=1.5;
    rightedge=1;   
    topedge=1;
    bottomedge=1.5;
end

if nargin < 5
    spacex=0.25;
    spacey=0.25;
end

plotheight=Ly*subplotsy + spacey*(subplotsy-1) + topedge + bottomedge;
plotwidth=Lx*subplotsx + spacex*(subplotsx-1) + leftedge + rightedge;

f=figure('visible','on');
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
sub_pos=ap_subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

hax = [];
for j = 1:subplotsx
    for k = 1:subplotsy
        ax = axes('position',sub_pos{j,subplotsy-k+1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
        hax = [hax; ax];
    end
end