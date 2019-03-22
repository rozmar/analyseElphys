function persistent_generatefigures_plot_APwaves(APwaves,valtozok)
if isfield(valtozok,'emphasize')
    axonalcolor=[1 .5 .5];
    somaticcolor=[.3 .3 .3];
    if isfield(valtozok.emphasize,'Window')
        ephidx= [APwaves.maxtime]>=valtozok.emphasize.Window(1)&[APwaves.maxtime]<=valtozok.emphasize.Window(2);
    else
        ephidx=false(size([APwaves.maxtime]));
    end
    if isfield(valtozok.emphasize,'idx')
        ephidx(valtozok.emphasize.idx)=true;
    end
    if ~isfield(valtozok.emphasize,'linewidth')
        valtozok.emphasize.linewidth=valtozok.linewidth;
    end 
else
    axonalcolor=[1 0 0];
    somaticcolor=[0 0 0];
end
if isield(valtozok,'dpi')
    dpi=valtozok.dpi;
else
    dpi=900;
end
if isield(valtozok,'dpi')
    dpi=valtozok.dpi;
else
    dpi=900;
end
if isfield(valtozok,'axeswidth')
    axesvastagsag=valtozok.axeswidth;
else
    axesvastagsag=1;
end
figure(1)
clf
hold on
needed=[APwaves.stimulated]&[APwaves.somaticAP];
plot([APwaves(needed).depolrate],[APwaves(needed).threshv]*1000,'ko','MarkerFaceColor',[0 0 0],'MarkerSize',valtozok.markersize,'LineWidth',valtozok.linewidth)
needed=[APwaves.stimulated]&[APwaves.axonalAP];
plot([APwaves(needed).depolrate],[APwaves(needed).threshv]*1000,'ro','MarkerFaceColor',[1 0 0],'MarkerSize',valtozok.markersize,'LineWidth',valtozok.linewidth)
needed=~[APwaves.stimulated]&[APwaves.axonalAP];
plot([APwaves(needed).depolrate],[APwaves(needed).threshv]*1000,'ro','MarkerFaceColor',[1 1 1],'MarkerSize',valtozok.markersize,'LineWidth',valtozok.linewidth)
needed=~[APwaves.stimulated]&[APwaves.somaticAP];
plot([APwaves(needed).depolrate],[APwaves(needed).threshv]*1000,'ko','MarkerFaceColor',[1 1 1],'MarkerSize',valtozok.markersize,'LineWidth',valtozok.linewidth)
xlabel('Rate of depolarization before AP (mV/ms)')
ylabel('Threshold (mV)')
if isfield(valtozok,'emphasize')
    needed=[APwaves.stimulated]&[APwaves.somaticAP] & ephidx;
    plot([APwaves(needed).depolrate],[APwaves(needed).threshv]*1000,'ro','MarkerFaceColor',[0 0 0],'MarkerSize',valtozok.markersize+1,'LineWidth',valtozok.emphasize.linewidth)
    needed=[APwaves.stimulated]&[APwaves.axonalAP]& ephidx;
    plot([APwaves(needed).depolrate],[APwaves(needed).threshv]*1000,'ro','MarkerFaceColor',[1 0 0],'MarkerSize',valtozok.markersize+1,'LineWidth',valtozok.emphasize.linewidth)
    needed=~[APwaves.stimulated]&[APwaves.axonalAP]& ephidx;
    plot([APwaves(needed).depolrate],[APwaves(needed).threshv]*1000,'ro','MarkerFaceColor',[1 1 1],'MarkerSize',valtozok.markersize+1,'LineWidth',valtozok.emphasize.linewidth)
    needed=~[APwaves.stimulated]&[APwaves.somaticAP]& ephidx;
    plot([APwaves(needed).depolrate],[APwaves(needed).threshv]*1000,'ko','MarkerFaceColor',[1 1 1],'MarkerSize',valtozok.markersize+1,'LineWidth',valtozok.emphasize.linewidth)
end

set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Fontname',betutipus,'Units','normalized','Position',[.25 .25 .5 .5])
set(gcf,'PaperUnits','centimeters','PaperPositionMode','manual','PaperSize',[xcm/.5 ycm/.5]+2,'PaperPosition',[2 2 xcm/.5 ycm/.5])
saveas(gcf,[dirs.figuresdir,expname,'_voltage.pdf'])
print(gcf,[dirs.figuresdir,expname,'_voltage.jpg'],'-djpeg',['-r',num2str(dpi)])

valtozok.linewidth=1;
figure(2)
clf
sis=unique([APwaves.si]);
for i=1:length(sis)
    
    subplot(2,2,1)
    hold on
    needed=find([APwaves.si]==sis(i) & [APwaves.stimulated]&[APwaves.somaticAP]) ;
    plot([APwaves(needed).t],[APwaves(needed).v],'-','Color',somaticcolor,'LineWidth',valtozok.linewidth);
    subplot(2,2,3)
    hold on
    plot([APwaves(needed).v],[APwaves(needed).dv],'-','Color',somaticcolor,'LineWidth',valtozok.linewidth);
    subplot(2,2,1)
    needed=find([APwaves.si]==sis(i) & [APwaves.stimulated]&[APwaves.axonalAP]) ;
    plot([APwaves(needed).t],[APwaves(needed).v],'-','Color',axonalcolor,'LineWidth',valtozok.linewidth);
    subplot(2,2,3)
    plot([APwaves(needed).v],[APwaves(needed).dv],'-','Color',axonalcolor,'LineWidth',valtozok.linewidth);
    subplot(2,2,2)
    hold on
    needed=find([APwaves.si]==sis(i) & ~[APwaves.stimulated]&[APwaves.somaticAP]) ;
    plot([APwaves(needed).t],[APwaves(needed).v],'-','Color',somaticcolor,'LineWidth',valtozok.linewidth);
    subplot(2,2,4)
    hold on
    plot([APwaves(needed).v],[APwaves(needed).dv],'-','Color',somaticcolor,'LineWidth',valtozok.linewidth);
    subplot(2,2,2)
    needed=find([APwaves.si]==sis(i) & ~[APwaves.stimulated]&[APwaves.axonalAP]) ;
    plot([APwaves(needed).t],[APwaves(needed).v],'-','Color',axonalcolor,'LineWidth',valtozok.linewidth);
    subplot(2,2,4)
    plot([APwaves(needed).v],[APwaves(needed).dv],'-','Color',axonalcolor,'LineWidth',valtozok.linewidth);
    
    if isfield(valtozok,'emphasize')
        subplot(2,2,1)
        hold on
        needed=find([APwaves.si]==sis(i) & [APwaves.stimulated]&[APwaves.somaticAP]& ephidx) ;
        plot([APwaves(needed).t],[APwaves(needed).v],'k-','LineWidth',valtozok.emphasize.linewidth);
        subplot(2,2,3)
        hold on
        plot([APwaves(needed).v],[APwaves(needed).dv],'k-','LineWidth',valtozok.emphasize.linewidth);
        subplot(2,2,1)
        needed=find([APwaves.si]==sis(i) & [APwaves.stimulated]&[APwaves.axonalAP]& ephidx) ;
        plot([APwaves(needed).t],[APwaves(needed).v],'r-','LineWidth',valtozok.emphasize.linewidth);
        subplot(2,2,3)
        plot([APwaves(needed).v],[APwaves(needed).dv],'r-','LineWidth',valtozok.emphasize.linewidth);
        subplot(2,2,2)
        hold on
        needed=find([APwaves.si]==sis(i) & ~[APwaves.stimulated]&[APwaves.somaticAP]& ephidx) ;
        plot([APwaves(needed).t],[APwaves(needed).v],'k-','LineWidth',valtozok.emphasize.linewidth);
        subplot(2,2,4)
        hold on
        plot([APwaves(needed).v],[APwaves(needed).dv],'k-','LineWidth',valtozok.emphasize.linewidth);
        subplot(2,2,2)
        needed=find([APwaves.si]==sis(i) & ~[APwaves.stimulated]&[APwaves.axonalAP] & ephidx);
        plot([APwaves(needed).t],[APwaves(needed).v],'r-','LineWidth',valtozok.emphasize.linewidth);
        subplot(2,2,4)
        plot([APwaves(needed).v],[APwaves(needed).dv],'r-','LineWidth',valtozok.emphasize.linewidth);
    end
    
    
    
end
subplot(2,2,1)
axis tight
ylims(1,:)=get(gca,'Ylim');
xlims(1,:)=get(gca,'Xlim');
subplot(2,2,2)
axis tight
ylims(2,:)=get(gca,'Ylim');
xlims(2,:)=get(gca,'Xlim');
subplot(2,2,3)
axis tight
ylims(3,:)=get(gca,'Ylim');
xlims(3,:)=get(gca,'Xlim');
subplot(2,2,4)
axis tight
ylims(4,:)=get(gca,'Ylim');
xlims(4,:)=get(gca,'Xlim');

subplot(2,2,1);
title('Current injection')
xlim([min(xlims(1:2,1)),max(xlims(1:2,2))])
ylim([min(ylims(1:2,1)),max(ylims(1:2,2))])
ylabel('Voltage (mV)')
xlabel('Time (ms)')
box off
subplot(2,2,2);
title('Spontaneous')
xlim([min(xlims(1:2,1)),max(xlims(1:2,2))])
ylim([min(ylims(1:2,1)),max(ylims(1:2,2))])
ylabel('Voltage (mV)')
xlabel('Time (ms)')
box off
subplot(2,2,3);
xlim([min(xlims(3:4,1)),max(xlims(3:4,2))])
ylim([min(ylims(3:4,1)),max(ylims(3:4,2))])
ylabel('dV/dt (mV/ms)')
xlabel('Voltage (mV)')
box off
subplot(2,2,4);
xlim([min(xlims(3:4,1)),max(xlims(3:4,2))])
ylim([min(ylims(3:4,1)),max(ylims(3:4,2))])
ylabel('dV/dt (mV/ms)')
xlabel('Voltage (mV)')
box off