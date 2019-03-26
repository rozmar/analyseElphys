function persistent_generatefigures_plot_APwaves(xlsidx,dirs,xlsdata,valtozok,expname,APwaves)%(APwaves,valtozok)

if isfield(valtozok,'xcm')
    xcm=valtozok.xcm;
else
    xcm=8;
end
if isfield(valtozok,'ycm')
    ycm=valtozok.ycm;
else
    ycm=8;
end

if isfield(valtozok,'highlight')
    axonalcolor=[1 .5 .5];
    somaticcolor=[.3 .3 .3];
    if isfield(valtozok.highlight,'Window')
        ephidx= [APwaves.maxtime]>=valtozok.highlight.Window(1)&[APwaves.maxtime]<=valtozok.highlight.Window(2);
    else
        ephidx=false(size([APwaves.maxtime]));
    end
    if isfield(valtozok.highlight,'idx')
        ephidx(valtozok.highlight.idx)=true;
    end
    if ~isfield(valtozok.highlight,'linewidth')
        valtozok.highlight.linewidth=valtozok.linewidth;
    end 
else
    axonalcolor=[1 0 0];
    somaticcolor=[0 0 0];
end
if isfield(valtozok,'dpi')
    dpi=valtozok.dpi;
else
    dpi=900;
end
if isfield(valtozok,'fontsize')
    betumeret=valtozok.fontsize;
else
    betumeret=14;
end
if isfield(valtozok,'fonttype')
    betutipus=valtozok.fonttype;
else
    betutipus='Arial';
end
if isfield(valtozok,'axeswidth')
    axesvastagsag=valtozok.axeswidth;
else
    axesvastagsag=1;
end

renderer='painters';

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
xlabel('Rate of depolarization (mV/ms)')
ylabel('Threshold (mV)')
if isfield(valtozok,'highlight')
    needed=[APwaves.stimulated]&[APwaves.somaticAP] & ephidx;
    plot([APwaves(needed).depolrate],[APwaves(needed).threshv]*1000,'ko','MarkerFaceColor',[0 0 0],'MarkerSize',valtozok.markersize+1,'LineWidth',valtozok.highlight.linewidth)
    needed=[APwaves.stimulated]&[APwaves.axonalAP]& ephidx;
    plot([APwaves(needed).depolrate],[APwaves(needed).threshv]*1000,'ro','MarkerFaceColor',[1 0 0],'MarkerSize',valtozok.markersize+1,'LineWidth',valtozok.highlight.linewidth)
    needed=~[APwaves.stimulated]&[APwaves.axonalAP]& ephidx;
    plot([APwaves(needed).depolrate],[APwaves(needed).threshv]*1000,'ro','MarkerFaceColor',[1 1 1],'MarkerSize',valtozok.markersize+1,'LineWidth',valtozok.highlight.linewidth)
    needed=~[APwaves.stimulated]&[APwaves.somaticAP]& ephidx;
    plot([APwaves(needed).depolrate],[APwaves(needed).threshv]*1000,'ko','MarkerFaceColor',[1 1 1],'MarkerSize',valtozok.markersize+1,'LineWidth',valtozok.highlight.linewidth)
end
axis tight
set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Fontname',betutipus,'Units','normalized','Position',[.25 .25 .5 .5])
set(gcf,'PaperUnits','centimeters','PaperPositionMode','manual','PaperSize',[xcm/.5 ycm/.5]+2,'PaperPosition',[2 2 xcm/.5 ycm/.5])
set(gcf, 'Renderer', renderer);
saveas(gcf,[dirs.figuresdir,expname,'_thresh_depolrate.pdf'])
print(gcf,[dirs.figuresdir,expname,'_thresh_depolrate.jpg'],'-djpeg',['-r',num2str(dpi)])

% valtozok.linewidth=1;
figure(2)
clf

figure(3)
clf

figure(4)
clf

figure(5)
clf


sis=unique([APwaves.si]);
for i=1:length(sis)
    
    figure(2)
    hold on
    needed=find([APwaves.si]==sis(i) & [APwaves.stimulated]&[APwaves.somaticAP]) ;
    plot([APwaves(needed).t],[APwaves(needed).v],'-','Color',somaticcolor,'LineWidth',valtozok.linewidth);
    figure(4)
    hold on
    plot([APwaves(needed).v],[APwaves(needed).dv],'-','Color',somaticcolor,'LineWidth',valtozok.linewidth);
    figure(2)
    needed=find([APwaves.si]==sis(i) & [APwaves.stimulated]&[APwaves.axonalAP]) ;
    plot([APwaves(needed).t],[APwaves(needed).v],'-','Color',axonalcolor,'LineWidth',valtozok.linewidth);
    figure(4)
    plot([APwaves(needed).v],[APwaves(needed).dv],'-','Color',axonalcolor,'LineWidth',valtozok.linewidth);
    figure(3)
    hold on
    needed=find([APwaves.si]==sis(i) & ~[APwaves.stimulated]&[APwaves.somaticAP]) ;
    plot([APwaves(needed).t],[APwaves(needed).v],'-','Color',somaticcolor,'LineWidth',valtozok.linewidth);
    figure(5)
    hold on
    plot([APwaves(needed).v],[APwaves(needed).dv],'-','Color',somaticcolor,'LineWidth',valtozok.linewidth);
    figure(3)
    needed=find([APwaves.si]==sis(i) & ~[APwaves.stimulated]&[APwaves.axonalAP]) ;
    plot([APwaves(needed).t],[APwaves(needed).v],'-','Color',axonalcolor,'LineWidth',valtozok.linewidth);
    figure(5)
    plot([APwaves(needed).v],[APwaves(needed).dv],'-','Color',axonalcolor,'LineWidth',valtozok.linewidth);
    
    if isfield(valtozok,'highlight')
        figure(2)
        hold on
        needed=find([APwaves.si]==sis(i) & [APwaves.stimulated]&[APwaves.somaticAP]& ephidx) ;
        plot([APwaves(needed).t],[APwaves(needed).v],'k-','LineWidth',valtozok.highlight.linewidth);
        figure(4)
        hold on
        plot([APwaves(needed).v],[APwaves(needed).dv],'k-','LineWidth',valtozok.highlight.linewidth);
        figure(2)
        needed=find([APwaves.si]==sis(i) & [APwaves.stimulated]&[APwaves.axonalAP]& ephidx) ;
        plot([APwaves(needed).t],[APwaves(needed).v],'r-','LineWidth',valtozok.highlight.linewidth);
        figure(4)
        plot([APwaves(needed).v],[APwaves(needed).dv],'r-','LineWidth',valtozok.highlight.linewidth);
        figure(3)
        hold on
        needed=find([APwaves.si]==sis(i) & ~[APwaves.stimulated]&[APwaves.somaticAP]& ephidx) ;
        plot([APwaves(needed).t],[APwaves(needed).v],'k-','LineWidth',valtozok.highlight.linewidth);
        figure(5)
        hold on
        plot([APwaves(needed).v],[APwaves(needed).dv],'k-','LineWidth',valtozok.highlight.linewidth);
        figure(3)
        needed=find([APwaves.si]==sis(i) & ~[APwaves.stimulated]&[APwaves.axonalAP] & ephidx);
        plot([APwaves(needed).t],[APwaves(needed).v],'r-','LineWidth',valtozok.highlight.linewidth);
        figure(5)
        plot([APwaves(needed).v],[APwaves(needed).dv],'r-','LineWidth',valtozok.highlight.linewidth);
    end
    
    
    
end
figure(2)
axis tight
ylims(1,:)=get(gca,'Ylim');
xlims(1,:)=get(gca,'Xlim');
figure(3)
axis tight
ylims(2,:)=get(gca,'Ylim');
xlims(2,:)=get(gca,'Xlim');
figure(4)
axis tight
ylims(3,:)=get(gca,'Ylim');
xlims(3,:)=get(gca,'Xlim');
figure(5)
axis tight
ylims(4,:)=get(gca,'Ylim');
xlims(4,:)=get(gca,'Xlim');

figure(2);
% title('Current injection')
xlim([min(xlims(1:2,1)),max(xlims(1:2,2))])
ylim([min(ylims(1:2,1)),max(ylims(1:2,2))])
ylabel('Voltage (mV)')
xlabel('Time (ms)')
box off
figure(3);
% title('Spontaneous')
xlim([min(xlims(1:2,1)),max(xlims(1:2,2))])
ylim([min(ylims(1:2,1)),max(ylims(1:2,2))])
ylabel('Voltage (mV)')
xlabel('Time (ms)')
box off
figure(4);
xlim([min(xlims(3:4,1)),max(xlims(3:4,2))])
ylim([min(ylims(3:4,1)),max(ylims(3:4,2))])
ylabel('dV/dt (mV/ms)')
xlabel('Voltage (mV)')
box off
figure(5);
xlim([min(xlims(3:4,1)),max(xlims(3:4,2))])
ylim([min(ylims(3:4,1)),max(ylims(3:4,2))])
ylabel('dV/dt (mV/ms)')
xlabel('Voltage (mV)')
box off

figure(2)
set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Fontname',betutipus,'Units','normalized','Position',[.25 .25 .5 .5])
set(gcf,'PaperUnits','centimeters','PaperPositionMode','manual','PaperSize',[xcm/.5 ycm/.5]+2,'PaperPosition',[2 2 xcm/.5 ycm/.5])
set(gcf, 'Renderer', renderer);
saveas(gcf,[dirs.figuresdir,expname,'_apwaves_stimulated_v_t.pdf'])
print(gcf,[dirs.figuresdir,expname,'_apwaves_stimulated_v_t.jpg'],'-djpeg',['-r',num2str(dpi)])

figure(3)
set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Fontname',betutipus,'Units','normalized','Position',[.25 .25 .5 .5])
set(gcf,'PaperUnits','centimeters','PaperPositionMode','manual','PaperSize',[xcm/.5 ycm/.5]+2,'PaperPosition',[2 2 xcm/.5 ycm/.5])
set(gcf, 'Renderer', renderer);
saveas(gcf,[dirs.figuresdir,expname,'_apwaves_spontaneous_v_t.pdf'])
print(gcf,[dirs.figuresdir,expname,'_apwaves_spontaneous_v_t.jpg'],'-djpeg',['-r',num2str(dpi)])

figure(4)
set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Fontname',betutipus,'Units','normalized','Position',[.25 .25 .5 .5])
set(gcf,'PaperUnits','centimeters','PaperPositionMode','manual','PaperSize',[xcm/.5 ycm/.5]+2,'PaperPosition',[2 2 xcm/.5 ycm/.5])
set(gcf, 'Renderer', renderer);
saveas(gcf,[dirs.figuresdir,expname,'_apwaves_stimulated_dv_v.pdf'])
print(gcf,[dirs.figuresdir,expname,'_apwaves_stimulated_dv_v.jpg'],'-djpeg',['-r',num2str(dpi)])

figure(5)
set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Fontname',betutipus,'Units','normalized','Position',[.25 .25 .5 .5])
set(gcf,'PaperUnits','centimeters','PaperPositionMode','manual','PaperSize',[xcm/.5 ycm/.5]+2,'PaperPosition',[2 2 xcm/.5 ycm/.5])
set(gcf, 'Renderer', renderer);
saveas(gcf,[dirs.figuresdir,expname,'_apwaves_spontaneous_dv_v.pdf'])
print(gcf,[dirs.figuresdir,expname,'_apwaves_spontaneous_dv_v.jpg'],'-djpeg',['-r',num2str(dpi)])