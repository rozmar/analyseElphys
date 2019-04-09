function persistent_generatefigures_plotselectedtimeinterval(xlsidx,dirs,xlsdata,valtozok,expname,additionaldata)

if isfield(valtozok,'xcm')
    xcm=valtozok.xcm;
else
    xcm=20;
end
if isfield(valtozok,'xcm_blowup')
    xcm_blowup=valtozok.xcm_blowup;
else
    xcm_blowup=xcm;
end
if isfield(valtozok,'ycm')
    ycm=valtozok.ycm;
else
    ycm=10;
end
if isfield(valtozok,'ycm_current')
    ycm_current=valtozok.ycm_current;
else
    ycm_current=ycm;
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
if isfield(valtozok,'voltagelinewidth')
    voltagelinewidth=valtozok.voltagelinewidth;
else
    voltagelinewidth=.5;
end
if isfield(valtozok,'currentlinewidth')
    currentlinewidth=valtozok.currentlinewidth;
else
    currentlinewidth=1;
end
if isfield(valtozok,'markersize')
    markersize=valtozok.markersize;
else
    markersize=1;
end
renderer='painters';
   %Usage: 
% ID='1608242rm_2_1,26_4'; 
% xlsidx=find(strcmp({xlsdata.ID},ID));
% valtozok=struct;
%valtozok.debugmode=0;
% valtozok.zerotime=76992;
% valtozok.xlimits=[0 60]+valtozok.zerotime;
% valtozok.xlimitstoexclude=[50,55;10,15]+valtozok.zerotime;
% valtozok.xlimitsblowup=[53, 59; 53 59] + valtozok.zerotime;
% valtozok.ylimitsblowup=[-80, 40;-80 -60];
% valtozok.ylimits=[-80 40;-30 50];
% expname=[ID,'  plateu current - voltage dependence'];
% valtozok.isiYlimits=[0 2];
% valtozok.freqYlimits=[.1 500];
% valtozok.freqYscale=[0.1, 1, 10, 100];
% valtozok.cutofffreq=5000;

% valtozok.threshold.threshold=10;
% valtozok.threshold.stepbacklength=.002;
% valtozok.threshold.aplength=.006;
% valtozok.threshold.displayonvoltagetrace=0;
% valtozok.threshold.displayonblowup=0;
% valtozok.threshold.movingstep=4;
% valtozok.threshold.threshmedianwindow=.15;
% valtozok.threshold.addbaselineval=1;

% valtozok.showevents=struct;
% valtozok.showevents(1).limits=[53,59]+valtozok.zerotime;
% valtozok.showevents(1).xlimits=[-.1 .5];
% valtozok.showevents(1).ylimits=[-80 40];
% valtozok.showevents(1).type='AP';
% valtozok.showevents(1).stimulated=true;

% valtozok.highlightaxonalspikes=1;
% valtozok.highlightaxonalspikes_timeback=.001;
% valtozok.highlightaxonalspikes_timeforward=.010;
% additionaldata=struct;
% additionaldata.eventdata=eventdata;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if size(valtozok.ylimits,1)>1
     valtozok.ylimitscurr=valtozok.ylimits(2,:);
     valtozok.ylimits=valtozok.ylimits(1,:);
 end


xinch=xcm/2.54;
yinch=ycm/2.54;
xsize=dpi*xinch;
ysize=dpi*yinch;


fname=[xlsdata(xlsidx).ID,'.mat'];
if nargin>5 & isfield(additionaldata,'eventdata')
    eventdata=additionaldata.eventdata;
else
    if ~isempty(dir(([dirs.eventdir,'sorted/',fname])));
        load([dirs.eventdir,'sorted/',fname],'eventdata');
    else
        eventdata=[];
    end
end
% eventdata=temp.eventdata;
load([dirs.bridgeddir,fname],'bridgeddata','stimdata');
% bridgeddata=temp.bridgeddata;
% stimdata=temp.stimdata;
if isfield(valtozok,'threshold')
    threshold=valtozok.threshold.threshold;
    stepbacklength=valtozok.threshold.stepbacklength;
    aplength=valtozok.threshold.aplength;
    stimstart=valtozok.xlimits(1);
    stimend=valtozok.xlimits(2);
    allapidxes=[eventdata.maxtime]>stimstart & [eventdata.maxtime]<stimend & strcmp({eventdata.type},'AP');%[stimapidxes,persistentapidxes];
    if isfield(valtozok,'xlimitstoexclude')
        for i=1:size(valtozok.xlimitstoexclude,1)
            allapidxes=allapidxes&~([eventdata.maxtime]>valtozok.xlimitstoexclude(i,1) & [eventdata.maxtime]<valtozok.xlimitstoexclude(i,2));
        end
    end
    allapidxes=find(allapidxes);
    for api=1:length(allapidxes)% calculating thresholds
        eventnum=allapidxes(api);
        sweepnum=eventdata(eventnum).sweepnum;
        ysweep=bridgeddata(sweepnum).y;
        ysweeporig=ysweep;
        
        si=eventdata(eventnum).si;
        starth=round(aplength/si);
        stepback=round(stepbacklength/si);
        
        y=NaN(1,starth*2+1);
        yorig=y;
        time=[-starth:starth]*si+eventdata(eventnum).maxtime+si;
        y(max((starth-eventdata(eventnum).maxh)+2,1):min(length(y),length(y)+length(ysweep)-(eventdata(eventnum).maxh+starth)))=ysweep(max(eventdata(eventnum).maxh-starth,1):min(eventdata(eventnum).maxh+starth,length(ysweep)));
        yorig(max((starth-eventdata(eventnum).maxh)+2,1):min(length(y),length(y)+length(ysweeporig)-(eventdata(eventnum).maxh+starth)))=ysweeporig(max(eventdata(eventnum).maxh-starth,1):min(eventdata(eventnum).maxh+starth,length(ysweeporig)));
        if valtozok.threshold.movingstep>1
            y=moving(y,valtozok.threshold.movingstep)'; % moving average!!
            yorig=moving(yorig,valtozok.threshold.movingstep)'; % moving average!!
        end
        stim(max((starth-eventdata(eventnum).maxh)+2,1):min(length(y),length(y)+length(ysweep)-(eventdata(eventnum).maxh+starth)))=stimdata(sweepnum).y(max(eventdata(eventnum).maxh-starth,1):min(eventdata(eventnum).maxh+starth,length(ysweep)));
        y(stim<stim(starth))=NaN;
        yorig(stim<stim(starth))=NaN;
        dy=diff(y)/si;
        y=mean([y(1:end-1);y(2:end)]);
        yorig=mean([yorig(1:end-1);yorig(2:end)]);
        time=mean([time(1:end-1);time(2:end)]);
        thresh=starth;
        while thresh>=1 & (y(thresh)>0 | y(thresh)+0.01>y(starth))
            thresh=thresh-1;
        end
%         while  thresh>=stepback+1 & nanmax(dy(thresh-stepback:thresh-1))>dy(thresh)
%             thresh=thresh-1;
%         end
        while thresh>=stepback+1 & nanmax(dy(thresh-stepback:thresh-1))>threshold
            thresh=thresh-1;
        end
        eventdata(eventnum).threshy=yorig(thresh);
        eventdata(eventnum).thresht=time(thresh);
    end
%      for api=1:length(allapidxes)
%          idx = abs([eventdata(allapidxes).thresht]-eventdata(allapidxes(api)).thresht)<valtozok.threshold.threshmedianwindow;
%            eventdata(allapidxes(api)).medianthreshy=median([eventdata(allapidxes(idx)).threshy]);
%      end
     threshtimestep=valtozok.threshold.threshmedianwindow/10;
     medianthreshtimevector=[valtozok.xlimits(1):threshtimestep:round(diff(valtozok.xlimits)/threshtimestep)*threshtimestep+valtozok.xlimits(1)]-valtozok.zerotime;
     medianthreshvalues=NaN(size(medianthreshtimevector));
     for i=1:length(medianthreshtimevector)
         idx = abs([eventdata(allapidxes).thresht]-valtozok.zerotime-medianthreshtimevector(i))<valtozok.threshold.threshmedianwindow;
         medianthreshvalues(i)=nanmedian([eventdata(allapidxes(idx)).threshy]);
     end
     todel=isnan(medianthreshvalues);
     medianthreshtimevector(todel)=[];
      medianthreshvalues(todel)=[];
%     channel=num2str(xlsdata(xlsnum).Channel);
%     sweeps=find([bridgeddata.realtime]>=stimstart & [bridgeddata.realtime]<=stimend & strcmp({bridgeddata.channellabel},['Vmon-',channel]));
%     sweepdata=struct;
    %          figure(1)
    %         clf
    %         hold on
%              for sweepi=1:length(sweeps)
%                  sweepnum=sweeps(sweepi);
%                  sweepdata(sweepi).y=bridgeddata(sweepnum).y;
%                  sweepdata(sweepi).x=[1:length(sweepdata(sweepi).y)]*bridgeddata(sweepnum).si+bridgeddata(sweepnum).realtime;
%                  plot([sweepdata(sweepi).x],[sweepdata(sweepi).y],'k-')
%              end
    %         plot([eventdata(allapidxes).thresht],[eventdata(allapidxes).threshy],'ro','LineWidth',2,'MarkerSize',6)
    %         pause
    
end

%%
neededidx=([bridgeddata.realtime]>valtozok.xlimits(1) &[bridgeddata.realtime]<valtozok.xlimits(2));
if isfield(valtozok,'xlimitstoexclude')
        for i=1:size(valtozok.xlimitstoexclude,1)
            neededidx=neededidx &~([bridgeddata.realtime]>valtozok.xlimitstoexclude(i,1) &[bridgeddata.realtime]<valtozok.xlimitstoexclude(i,2));
        end
end
neededidx=find(neededidx);
if neededidx(1)~=1
    neededidx=[neededidx(1)-1,neededidx];
end
bridgedsweeps=struct;
APsnow=struct;
for sweep=1:length(neededidx)
    if valtozok.cutofffreq(1)>0
        if length(valtozok.cutofffreq)==1
            [b,a]=butter(1,valtozok.cutofffreq/(1/bridgeddata(neededidx(sweep)).si)/2,'low');
        else
            [b,a]=butter(1,valtozok.cutofffreq/(1/bridgeddata(neededidx(sweep)).si)/2,'bandpass');
        end
    end
    bridgedsweeps(sweep).starttime=bridgeddata(neededidx(sweep)).realtime;
    if valtozok.cutofffreq>0
        bridgedsweeps(sweep).y=filtfilt(b,a,bridgeddata(neededidx(sweep)).y)*1000;
    else
        bridgedsweeps(sweep).y=bridgeddata(neededidx(sweep)).y*1000;
    end
    bridgedsweeps(sweep).stim=stimdata(neededidx(sweep)).y*10^12;
    bridgedsweeps(sweep).x=[bridgeddata(neededidx(sweep)).realtime:bridgeddata(neededidx(sweep)).si:bridgeddata(neededidx(sweep)).realtime+(length(bridgedsweeps(sweep).y)-1)*bridgeddata(neededidx(sweep)).si]-valtozok.zerotime;
    bridgedsweeps(sweep).sweepnum=neededidx(sweep);
    bridgedsweeps(sweep).RS=stimdata(neededidx(sweep)).RS*1e-6;
    if isfield(valtozok,'threshold') & isfield(valtozok.threshold,'addbaselineval') & valtozok.threshold.addbaselineval==1 % baseline calculation
        idx=find(bridgedsweeps(sweep).stim~=bridgedsweeps(sweep).stim(1),1,'first');
        if isempty(idx)
            idx=length(bridgedsweeps(sweep).y);
        end
        bridgedsweeps(sweep).baselineval=nanmin( bridgedsweeps(sweep).y(1:idx));
        bridgedsweeps(sweep).baselinetime=bridgedsweeps(sweep).x(idx);
    end
end
figure(1)
clf
hold on
for sweep=1:length(bridgedsweeps)
    plot(bridgedsweeps(sweep).x,bridgedsweeps(sweep).y,['k-'],'LineWidth',voltagelinewidth)
    
    %%
    if ~isfield(valtozok,'highlightaxonalspikes') | valtozok.highlightaxonalspikes==1
%         eventdatanow=eventdata([strcmp({eventdata.type},'AP')]& ~[eventdata.stimulated] & [eventdata.sweepnum]==bridgedsweeps(sweep).sweepnum);
            eventdatanow=eventdata([eventdata.axonalAP]==1& [eventdata.sweepnum]==bridgedsweeps(sweep).sweepnum & ~[eventdata.stimulated]) ;
            si=mode(diff(bridgedsweeps(sweep).x));
            stepback=round(valtozok.highlightaxonalspikes_timeback/si);
            stepforward=round(valtozok.highlightaxonalspikes_timeforward/si);
        for i=1:length(eventdatanow)
%             idxes=eventdatanow(i).onseth:eventdatanow(i).endh;
            idxes=[-stepback:stepforward]+eventdatanow(i).maxh;
            idxes(idxes<=0 | idxes>length(bridgedsweeps(sweep).y))=[];
            plot(bridgedsweeps(sweep).x(idxes),bridgedsweeps(sweep).y(idxes),['r-'],'LineWidth',voltagelinewidth)
        end
    end
end
if isfield(valtozok,'threshold') & valtozok.threshold.displayonvoltagetrace==1
    plot([eventdata(allapidxes).thresht]-valtozok.zerotime,[eventdata(allapidxes).threshy]*1000,'ro','LineWidth',1,'MarkerSize',5)
    %       plot([eventdata(allapidxes).thresht]-valtozok.zerotime,[eventdata(allapidxes).medianthreshy]*1000,'rs-','LineWidth',2,'MarkerSize',5,'MarkerFaceColor',[1 0 0])
    plot(medianthreshtimevector,medianthreshvalues*1000,'r-','LineWidth',2,'MarkerSize',5,'MarkerFaceColor',[1 0 0])
    if isfield(valtozok,'threshold') & isfield(valtozok.threshold,'addbaselineval') & valtozok.threshold.addbaselineval==1
        plot([bridgedsweeps.baselinetime],[bridgedsweeps.baselineval],'bo-','LineWidth',2,'MarkerSize',5,'MarkerFaceColor',[0 0 1])
    end
end
%%
if isfield(valtozok,'drugwashin') & ~isempty(xlsdata(xlsidx).drugdata)
    drugtimes=[[xlsdata(xlsidx).drugdata.DrugWashinTime];[xlsdata(xlsidx).drugdata.DrugWashoutTime]]-valtozok.zerotime;
    uniquedrugtimes=unique(drugtimes','rows');
    drugwashinline_ystart=valtozok.drugwashin.drugwashinline_ystart;
    drugwashinline_ystep=valtozok.drugwashin.drugwashinline_ystep;
    drugwashinline_linewidth=valtozok.drugwashin.drugwashinline_linewidth;
    for drugi=1:size(uniquedrugtimes,1)
        yvalnow=drugwashinline_ystart+drugwashinline_ystep*(drugi-1);
        if uniquedrugtimes(drugi,1)<=valtozok.xlimits(2)-valtozok.zerotime
            plot([max(uniquedrugtimes(drugi,1),valtozok.xlimits(1)-valtozok.zerotime),min(uniquedrugtimes(drugi,2),valtozok.xlimits(2)-valtozok.zerotime)],[yvalnow,yvalnow],'k-','LineWidth',drugwashinline_linewidth);
        end
    end
end
%%

axis tight
xlim(valtozok.xlimits-valtozok.zerotime)
if diff(valtozok.ylimits)>0
    ylim(valtozok.ylimits)
end
valtozok.ylimits=get(gca,'Ylim');
% %%
% set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Fontname',betutipus,'Position',[1/xcm 1/ycm 1-2/xcm 1-2/ycm])
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 xsize/dpi ysize/dpi])
% print(gcf,[dirs.figuresdir,expname,'_voltage.jpg'],'-djpeg',['-r',num2str(dpi)])
% saveas(gcf,[dirs.figuresdir,expname,'_voltage.pdf'])
%%
xlabel('Time (s)')
ylabel('Voltage (mV)')

%%
if isfield(valtozok,'axis')
    h = gca;
    if valtozok.axis.voltage_x
        h.XAxis.Visible = 'on';
    else
        h.XAxis.Visible = 'off';
    end
    if valtozok.axis.voltage_y
        h.YAxis.Visible = 'on';
    else
        h.YAxis.Visible = 'of';
    end
end
%%
set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Fontname',betutipus,'Units','normalized','Position',[.25 .25 .5 .5])
set(findobj(gcf,'type','text'),'fontsize',betumeret,'Fontname',betutipus)
set(gcf,'PaperUnits','centimeters','PaperPositionMode','manual','PaperSize',[xcm/.5 ycm/.5]+2,'PaperPosition',[2 2 xcm/.5 ycm/.5])
set(gcf, 'Renderer', renderer);
% saveas(gcf,[dirs.figuresdir,expname,'_voltage.pdf'])
print(gcf,[dirs.figuresdir,expname,'_voltage.jpg'],'-djpeg',['-r',num2str(dpi)])

%%


if isfield(valtozok,'debugmode')&valtozok.debugmode==1
    return
end
if sum(diff(valtozok.xlimitsblowup(1,:)))>0
    if isfield(valtozok,'axis') & isfield(valtozok.axis,'voltage_blowup_x')
        h = gca;
        if valtozok.axis.voltage_blowup_x
            h.XAxis.Visible = 'on';
        else
            h.XAxis.Visible = 'off';
        end
        if valtozok.axis.voltage_blowup_y
            h.YAxis.Visible = 'on';
        else
            h.YAxis.Visible = 'of';
        end
    end
    if isfield(valtozok,'threshold') & valtozok.threshold.displayonblowup==1
      plot([eventdata(allapidxes).thresht]-valtozok.zerotime,[eventdata(allapidxes).threshy]*1000,'ro','LineWidth',1,'MarkerSize',5)
      plot(medianthreshtimevector,medianthreshvalues*1000,'r-','LineWidth',2,'MarkerSize',5,'MarkerFaceColor',[1 0 0])
       if isfield(valtozok,'threshold') & isfield(valtozok.threshold,'addbaselineval') & valtozok.threshold.addbaselineval==1
           plot([bridgedsweeps.baselinetime],[bridgedsweeps.baselineval],'bo-','LineWidth',2,'MarkerSize',5,'MarkerFaceColor',[0 0 1])
       end
    end
    for blowupi=1:size(valtozok.xlimitsblowup,1)
        
        xlim(valtozok.xlimitsblowup(blowupi,:)-valtozok.zerotime);
        if diff(valtozok.ylimitsblowup(blowupi,:))>0
            ylim(valtozok.ylimitsblowup(blowupi,:))
        else
            ylim(valtozok.ylimits)
        end
        
        set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Fontname',betutipus,'Units','normalized','Position',[.25 .25 .5 .5])
        set(findobj(gcf,'type','text'),'fontsize',betumeret,'Fontname',betutipus)
        set(gcf,'PaperUnits','centimeters','PaperPositionMode','manual','PaperSize',[xcm_blowup/.5 ycm/.5]+2,'PaperPosition',[2 2 xcm_blowup/.5 ycm/.5])
        set(gcf, 'Renderer', renderer);
        saveas(gcf,[dirs.figuresdir,expname,'_voltage_blowup_',num2str(blowupi),'_.pdf'])
        print(gcf,[dirs.figuresdir,expname,'_voltage_blowup_',num2str(blowupi),'_.jpg'],'-djpeg',['-r',num2str(dpi)])
    end
end


ycmnow=ycm_current;%/2;
xinch=xcm/2.54;
yinch=ycmnow/2.54;
xsize=dpi*xinch;
ysize=dpi*yinch;
figure(1)
clf
hold on
for sweep=1:length(bridgedsweeps)
    plot(bridgedsweeps(sweep).x,bridgedsweeps(sweep).stim,'k-','LineWidth',currentlinewidth)
end
axis tight
xlim(valtozok.xlimits-valtozok.zerotime)
if isfield(valtozok,'ylimitscurr');
    ylim(valtozok.ylimitscurr)
end
% set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Fontname',betutipus,'Position',[1/xcm 1/ycmnow 1-2/xcm 1-2/ycmnow])
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 xsize/dpi ysize/dpi])
% print(gcf,[dirs.figuresdir,expname,'_current.jpg'],'-djpeg',['-r',num2str(dpi)])
% saveas(gcf,[dirs.figuresdir,expname,'_current.pdf'])
xlabel('Time (s)')
ylabel('Injected current (pA)')
if isfield(valtozok,'axis')
    h = gca;
    if valtozok.axis.current_x
        h.XAxis.Visible = 'on';
    else
        h.XAxis.Visible = 'off';
    end
    if valtozok.axis.current_y
        h.YAxis.Visible = 'on';
    else
        h.YAxis.Visible = 'of';
    end
end
set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Fontname',betutipus,'Units','normalized','Position',[.25 .25 .5 .5])
set(findobj(gcf,'type','text'),'fontsize',betumeret,'Fontname',betutipus)
set(gcf,'PaperUnits','centimeters','PaperPositionMode','manual','PaperSize',[xcm/.5 ycm_current/.5]+2,'PaperPosition',[2 2 xcm/.5 ycm_current/.5])
set(gcf, 'Renderer', renderer);
saveas(gcf,[dirs.figuresdir,expname,'_current.pdf'])
print(gcf,[dirs.figuresdir,expname,'_current.jpg'],'-djpeg',['-r',num2str(dpi)])
if ~isempty(eventdata)
eventdataold=eventdata;
eventdata=eventdata([strcmp({eventdata.type},'AP')]);
eventdata=eventdata([eventdata.maxtime]>=valtozok.xlimits(1) & [eventdata.maxtime]<=valtozok.xlimits(2));
%%
APpers=eventdata([eventdata.stimulated]==0);
APstim=eventdata([eventdata.stimulated]==1);

isis=[0,diff([APpers.maxtime])];
isisstim=[0,diff([APstim.maxtime])];

figure(2)
clf
hold on

if length(isis)>1
    plot([APpers.maxtime]-valtozok.zerotime,isis,'ro','MarkerSize',markersize,'MarkerFaceColor',[1 0 0])
end
if length(isisstim)>1
    plot([APstim.maxtime]-valtozok.zerotime,isisstim,'ko','MarkerSize',markersize,'MarkerFaceColor',[0 0 0])
end
axis tight
if diff(valtozok.isiYlimits)>0
    ylim(valtozok.isiYlimits)
end
xlim(valtozok.xlimits-valtozok.zerotime)
box off
% set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Fontname',betutipus,'Position',[1/xcm 1/ycm 1-2/xcm 1-2/ycm])
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 xsize/dpi ysize/dpi])
if isfield(valtozok,'axis')
    h = gca;
    if valtozok.axis.freq_x
        h.XAxis.Visible = 'on';
    else
        h.XAxis.Visible = 'off';
    end
    if valtozok.axis.freq_y
        h.YAxis.Visible = 'on';
    else
        h.YAxis.Visible = 'of';
    end
end
set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Fontname',betutipus,'Units','normalized','Position',[.25 .25 .5 .5])
set(findobj(gcf,'type','text'),'fontsize',betumeret,'Fontname',betutipus)
set(gcf,'PaperUnits','centimeters','PaperPositionMode','manual','PaperSize',[xcm/.5 ycm/.5]+2,'PaperPosition',[2 2 xcm/.5 ycm/.5])
set(gcf, 'Renderer', renderer);
print(gcf,[dirs.figuresdir,expname,'_ISI.jpg'],'-djpeg',['-r',num2str(dpi)])
% betumeret=14;
% axesvastagsag=1;
figure(4)
clf
if length(isis)>1
    semilogy([APpers.maxtime]-valtozok.zerotime,1./isis,'ro','MarkerSize',markersize,'MarkerFaceColor',[1 0 0])
end
hold on
if length(isisstim)>1
    semilogy([APstim.maxtime]-valtozok.zerotime,1./isisstim,'ko','MarkerSize',markersize,'MarkerFaceColor',[0 0 0])
end
% xlabel('time(s)')
% ylabel('Spont. Frequency (Hz)')
axis tight
xlim(valtozok.xlimits-valtozok.zerotime)
if diff(valtozok.freqYlimits)>0
    ylim(valtozok.freqYlimits)
end
if sum(valtozok.freqYscale)>0
    set(gca,'Ytick',valtozok.freqYscale)
end
box off

if isfield(valtozok,'axis')
    h = gca;
    if valtozok.axis.freq_x
        h.XAxis.Visible = 'on';
    else
        h.XAxis.Visible = 'off';
    end
    if valtozok.axis.freq_y
        h.YAxis.Visible = 'on';
    else
        h.YAxis.Visible = 'of';
    end
end

% set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Fontname',betutipus,'Position',[1/xcm 1/ycmnow 1-2/xcm 1-2/ycmnow])
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 xsize/dpi ysize/dpi])
set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Fontname',betutipus,'Units','normalized','Position',[.25 .25 .5 .5])
set(findobj(gcf,'type','text'),'fontsize',betumeret,'Fontname',betutipus)
set(gcf,'PaperUnits','centimeters','PaperPositionMode','manual','PaperSize',[xcm/.5 ycm/.5]+2,'PaperPosition',[2 2 xcm/.5 ycm/.5])
set(gcf, 'Renderer', renderer);
print(gcf,[dirs.figuresdir,expname,'_freq.jpg'],'-djpeg',['-r',num2str(dpi)])
%%
eventdata=eventdataold;
figure(5)
clf
hold on
if isfield(valtozok,'threshold')
    persap=allapidxes(~[eventdata(allapidxes).stimulated]);
    stimap=allapidxes([eventdata(allapidxes).stimulated]);
    plot([eventdata(stimap).thresht]-valtozok.zerotime,[eventdata(stimap).threshy]*1000,'ko','MarkerSize',3,'MarkerFaceColor',[0 0 0])
    plot([eventdata(persap).thresht]-valtozok.zerotime,[eventdata(persap).threshy]*1000,'ro','MarkerSize',3,'MarkerFaceColor',[1 0 0])
    
%     plot(medianthreshtimevector,medianthreshvalues*1000,'r-','LineWidth',2,'MarkerSize',5,'MarkerFaceColor',[1 0 0])
    if isfield(valtozok,'threshold') & isfield(valtozok.threshold,'addbaselineval') & valtozok.threshold.addbaselineval==1
        plot([bridgedsweeps.baselinetime],[bridgedsweeps.baselineval],'bo-','LineWidth',2,'MarkerSize',5,'MarkerFaceColor',[0 0 1])
    end
    axis tight
    xlim(valtozok.xlimits-valtozok.zerotime);

    box off
    set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Fontname',betutipus,'Position',[1/xcm 1/ycmnow 1-2/xcm 1-2/ycmnow])
    set(findobj(gcf,'type','text'),'fontsize',betumeret,'Fontname',betutipus)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 xsize/dpi ysize/dpi])
    print(gcf,[dirs.figuresdir,expname,'_thresh.jpg'],'-djpeg',['-r',num2str(dpi)])
end
end
%%
% betumeret=14;
% betutipus='Arial';
% axesvastagsag=1;

xinch=xcm/2.54;
yinch=ycm/2.54;
xsize=dpi*xinch;
ysize=dpi*yinch;

if isfield(valtozok,'showevents')
    for i=1:length(valtozok.showevents)
        neededevents=find([eventdata.maxtime]>valtozok.showevents(i).limits(1) & [eventdata.maxtime]<valtozok.showevents(i).limits(2) & strcmp({eventdata.type},valtozok.showevents(i).type) & [eventdata.stimulated]==valtozok.showevents(i).stimulated);
        figure(313442)
        clf
        hold on
        for eventi=1:length(neededevents)
            sweepnum=eventdata(neededevents(eventi)).sweepnum;
            maxh=eventdata(neededevents(eventi)).maxh;
            si=eventdata(neededevents(eventi)).si;
            if valtozok.cutofffreq>0
                [b,a]=butter(3,valtozok.cutofffreq/(1/si)/2,'low');
                y=filtfilt(b,a,bridgeddata(sweepnum).y)*1000;
            else
                y=bridgeddata(sweepnum).y*1000;
            end
            
            
            time=[1:length(y)]*si-maxh*si;
            plot(time,y,'k-','LineWidth',1)
        end
        xlim(valtozok.showevents(i).xlimits);
        if diff(valtozok.showevents(i).ylimits)>0
            ylim(valtozok.showevents(i).ylimits)
        end
        set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Fontname',betutipus,'Position',[1/xcm 1/ycm 1-2/xcm 1-2/ycm])
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 xsize/dpi ysize/dpi])
    set(gcf, 'Renderer', renderer);
    print(gcf,[dirs.figuresdir,expname,'_events_',num2str(i),'.jpg'],'-djpeg',['-r',num2str(dpi)])
    end
    
end