%% recording statistics
anaesthgroups=unique({xlsdata.anaesthesia},'stable');
baselineSDbinsize=5;
baselineSDbinnum=round(30*60/baselineSDbinsize);
recstats=struct;
files={xlsdata.ID};
for filei=1:length(files)
    if xlsdata(filei).field==0
        if isempty(fieldnames(recstats))
            NEXT=1;
        else
            NEXT=length(recstats)+1;
        end
        ID=files{filei};
        load([dirs.bridgeddir,ID],'bridgeddata','lightdata','stimdata');
        recstats(NEXT).ID=ID;
        recstats(NEXT).RS=nanmedian([lightdata.RS]);
        recstats(NEXT).anaesth=xlsdata(filei).anaesthesia;
        if lightdata(end).realtime>lightdata(1).realtime;
            recstats(NEXT).recordinglength=lightdata(end).realtime-lightdata(1).realtime;
        else
            recstats(NEXT).recordinglength=lightdata(end).realtime-lightdata(1).realtime+24*3600;
        end
        recstats(NEXT).recordinglength= recstats(NEXT).recordinglength+bridgeddata(end).si*length(bridgeddata(end).y);
        %baselineSD statistics
        baselineSD=NaN(baselineSDbinnum,1);
        if recstats(NEXT).recordinglength>=20*60
            for bini=1:baselineSDbinnum
                startt=(bini-1)*baselineSDbinsize+lightdata(1).realtime;
                endt=(bini)*baselineSDbinsize+lightdata(1).realtime;
                neededsweepnum=find([lightdata.realtime]>=startt & [lightdata.realtime]<=endt);
                sweepSD=NaN;
                if ~isempty(neededsweepnum)
                    if neededsweepnum(1)>1
                        neededsweepnum=[neededsweepnum(1)-1,neededsweepnum];
                    end
                    sweepSD=nan(size(neededsweepnum));
                    for sweepi=1:length(neededsweepnum)
                        sweepnum=neededsweepnum(sweepi);
                        if strcmp(stimdata(sweepnum).Amplifiermode,'C-Clamp') &  any(strfind(bridgeddata(sweepnum).channellabel,'Vmon')) & std(stimdata(sweepnum).y)==0
                            [b,a]=butter(1,500/(1/bridgeddata(sweepnum).si)/2,'low');
                            y=bridgeddata(sweepnum).y;
                            y=filtfilt(b,a,y);
                            x=[1:length(bridgeddata(sweepnum).y)]*bridgeddata(sweepnum).si+bridgeddata(sweepnum).realtime-bridgeddata(sweepnum).si;
                            neededidx=find(x>=startt & x<=endt);
                            sweepSD(sweepi)=nanstd(y(neededidx));
                        end
                    end
                end
                baselineSD(bini)=nanmedian(sweepSD);
            end
            %             figure(2)
            %             clf
            %             plot([1:baselineSDbinnum]*baselineSDbinsize/60,baselineSD,'ko-')
            %             pause
        end
        recstats(NEXT).baselineSD=baselineSD;
    end
end

baselineSDtimevector=[1:baselineSDbinnum]*baselineSDbinsize/60;
newbaselineSDtimevector=[1:1:30];
newbaselineSDtimevectorcenters=newbaselineSDtimevector-mode(diff(newbaselineSDtimevector));
newbaselineSDtimevectorsteps=mode(diff(newbaselineSDtimevector));
for i=1:length(recstats)
    recstats(i).newbaselineSD=[];
    for stepi=1:length(newbaselineSDtimevectorcenters)
        idx=find(baselineSDtimevector>=newbaselineSDtimevectorcenters(stepi)-newbaselineSDtimevectorsteps/2 & baselineSDtimevector<=newbaselineSDtimevectorcenters(stepi)+newbaselineSDtimevectorsteps/2);
        recstats(i).newbaselineSD(stepi)=nanmean(recstats(i).baselineSD(idx));
    end
    recstats(i).newbaselineSD=(recstats(i).newbaselineSD/nanmedian(recstats(i).newbaselineSD(find(newbaselineSDtimevector>15))))';
end

% time of the recordings
for i=1:length(recstats)
    recdate=recstats(i).ID(1:6);
    recstats(i).datenumber=datenum(str2num((['20',recdate(1:2)])),str2num(recdate(3:4)),str2num(recdate(5:6)));
     recstats(i).datestring=datestr(datenum(str2num((['20',recdate(1:2)])),str2num(recdate(3:4)),str2num(recdate(5:6))));
     datevector=datevec(datenum(str2num((['20',recdate(1:2)])),str2num(recdate(3:4)),str2num(recdate(5:6))));
     recstats(i).year=datevector(1);
     recstats(i).month=datevector(2);
     recstats(i).day=datevector(3);
end

%%
startdate=datenum(2017,01,01);
enddate=datenum(2019,01,01);
figure(1)
clf
subplot(3,2,1)
[nall,xbinn]=hist([recstats.recordinglength]/60,[2.5:5:62.5]);
[RSall,xbinRS]=hist([recstats.RS]/10^6,[5:10:95]);
[dateall,xbindate]=hist([recstats.datenumber],[startdate:7:enddate]);
xbinZ=[5:10:155];
xbinage=[0:10:366];
histstruct=struct;
for groupi=1:length(anaesthgroups)
    [histstruct(groupi).n,~]=hist([recstats(strcmp({recstats.anaesth},anaesthgroups{groupi})).recordinglength]/60,xbinn);
    histstruct(groupi).n=histstruct(groupi).n';
    [histstruct(groupi).RS,~]=hist([recstats(strcmp({recstats.anaesth},anaesthgroups{groupi})).RS]/10^6,xbinRS);
    histstruct(groupi).RS=histstruct(groupi).RS';
    histstruct(groupi).Z=[];
    histstruct(groupi).age=[];
    [histstruct(groupi).date,~]=hist([recstats(strcmp({recstats.anaesth},anaesthgroups{groupi})).datenumber],xbindate);
    histstruct(groupi).date=histstruct(groupi).date';
    
    for xlsi=1:length(xlsdata)
        if isnumeric( xlsdata(xlsi).locationz) & strcmp(xlsdata(xlsi).anaesthesia,anaesthgroups{groupi}) & xlsdata(xlsi).field==0
            histstruct(groupi).Z=[histstruct(groupi).Z, xlsdata(xlsi).locationz];
            histstruct(groupi).age=[histstruct(groupi).age, xlsdata(xlsi).age];
            
        end
    end
    [histstruct(groupi).Zhist,~]=hist(histstruct(groupi).Z,xbinZ);
    histstruct(groupi).Zhist=histstruct(groupi).Zhist';
    [histstruct(groupi).Agehist,~]=hist(histstruct(groupi).age,xbinage);
    histstruct(groupi).Agehist=histstruct(groupi).Agehist';
end
[hAxes,hBar,hLine] = plotyy(xbinn,[histstruct.n],xbinn,cumsum([histstruct.n]),'bar','plot');

set(hBar,'BarLayout','stacked')
% set(gca,'box','off')
ylim([0,max(sum([histstruct.n]'))]);
% set(gca,'Ytick',[0:max(sum([histstruct.n]'))]);
set(hAxes(1), 'YLimMode', 'Auto');
set(hAxes(1),'Ytickmode','auto');
set(hAxes(2),'Ytickmode','auto');
set(hAxes(2), 'YLimMode', 'Auto');
set(hAxes(2),'XTick',[],'XTickLabel',[]);
set(hAxes(1),'Xlim',[0 max(xbinn)]);
set(hAxes(2),'Xlim',[0 max(xbinn)]);
% axes(hAxes(2))
% axis tight
colororder=get(hAxes,'colororder');
for groupi=2:length(anaesthgroups)
    set(hBar(groupi),'FaceColor',colororder{1}(groupi,:));
end
set(hLine,'LineWidth',2)

set(gca,'Xtick',[0:5:60])



xlabel('recording length (min)')
legend(anaesthgroups)
% axes(hAxes(1))
ylabel('# of cells')
% title('recording length')


subplot(3,2,2)
[hAxes,hBar,hLine] = plotyy(xbinRS,[histstruct.RS],xbinRS,cumsum([histstruct.RS]),'bar','plot');
set(hBar,'BarLayout','stacked')
set(gca,'box','off')
% ylim([0,max(sum([histstruct.RS]'))]);
set(hAxes(1), 'YLimMode', 'Auto');
set(hAxes(1),'Ytickmode','auto');
set(hAxes(2),'Ytickmode','auto');
set(hAxes(2), 'YLimMode', 'Auto');
set(hAxes(2),'XTick',[],'XTickLabel',[]);
set(hAxes(1),'Xlim',[0 max(xbinRS)]);
set(hAxes(2),'Xlim',[0 max(xbinRS)]);
colororder=get(hAxes,'colororder');
for groupi=2:length(anaesthgroups)
    set(hBar(groupi),'FaceColor',colororder{1}(groupi,:));
end
set(hLine,'LineWidth',2)
% legend(anaesthgroups)
set(gca,'Xtick',[0:10:90])
xlabel('median RS (MOhm)')
ylabel('# of cells')
% title('median RS')


subplot(3,2,3)
[hAxes,hBar,hLine] = plotyy(xbinZ,[histstruct.Zhist],xbinZ,cumsum([histstruct.Zhist]),'bar','plot');
set(hBar,'BarLayout','stacked')
set(gca,'box','off')
% ylim([0,max(sum([histstruct.RS]'))]);
set(hAxes(1), 'YLimMode', 'Auto');
set(hAxes(1),'Ytickmode','auto');
set(hAxes(2),'Ytickmode','auto');
set(hAxes(2), 'YLimMode', 'Auto');
set(hAxes(2),'XTick',[],'XTickLabel',[]);
set(hAxes(1),'Xlim',[0 max(xbinZ)]);
set(hAxes(2),'Xlim',[0 max(xbinZ)]);
colororder=get(hAxes,'colororder');
for groupi=2:length(anaesthgroups)
    set(hBar(groupi),'FaceColor',colororder{1}(groupi,:));
end
set(hLine,'LineWidth',2)
% legend(anaesthgroups)
set(gca,'Xtick',[0:10:140])
xlabel('Depth from pia (microns)')
ylabel('# of cells')
% title('Depth from the surface of the brain')

subplot(3,2,4)
cla
[hAxes,hBar,hLine] = plotyy(xbindate,[histstruct.date],xbindate,cumsum([histstruct.date]),'bar','plot');
set(hBar,'BarLayout','stacked')
set(gca,'box','off')
axis tight
% ylim([0,max(sum([histstruct.RS]'))]);
set(hAxes(1), 'YLimMode', 'Auto');
set(hAxes(1),'Ytickmode','auto');
set(hAxes(2),'Ytickmode','auto');
set(hAxes(2), 'YLimMode', 'Auto');
set(hAxes(2),'XTick',[],'XTickLabel',[]);
set(hAxes(1),'Xlim',[min(xbindate) max(xbindate)]);
set(hAxes(2),'Xlim',[min(xbindate) max(xbindate)]);
colororder=get(hAxes,'colororder');
for groupi=2:length(anaesthgroups)
    set(hBar(groupi),'EdgeColor',[0 0 0],'FaceColor',colororder{1}(groupi,:));
end
set(hLine,'LineWidth',2)
% legend(anaesthgroups)
% set(gca,'Xtick',[0:10:140])
xlabel('Week of recording')
ylabel('# of cells')
newlabel=datestr(cellfun(@str2num, get(gca,'XtickLabels'))*10^5);
set(gca,'XtickLabels',newlabel)
% disp(datestr(cellfun(@str2num, get(gca,'XtickLabels'))))
% title('Depth from the surface of the brain')

subplot(3,2,5)
cla
[hAxes,hBar,hLine] = plotyy(xbinage,[histstruct.Agehist],xbinage,cumsum([histstruct.Agehist]),'bar','plot');
set(hBar,'BarLayout','stacked')
set(gca,'box','off')
axis tight
% ylim([0,max(sum([histstruct.RS]'))]);
set(hAxes(1), 'YLimMode', 'Auto');
set(hAxes(1),'Ytickmode','auto');
set(hAxes(1),'Xlim',[min(xbinage) max(xbinage)]);
set(hAxes(2),'Xlim',[min(xbinage) max(xbinage)]);
set(hAxes(2),'Ytickmode','auto');
set(hAxes(2), 'YLimMode', 'Auto');
set(hAxes(2),'XTick',[],'XTickLabel',[]);
colororder=get(hAxes,'colororder');
for groupi=2:length(anaesthgroups)
    set(hBar(groupi),'EdgeColor',[0 0 0],'FaceColor',colororder{1}(groupi,:));
end
set(hLine,'LineWidth',2)
% legend(anaesthgroups)
% set(gca,'Xtick',[0:10:140])
xlabel('Age of animal (days)')
ylabel('# of cells')
% disp(datestr(cellfun(@str2num, get(gca,'XtickLabels'))))
% title('Depth from the surface of the brain')


subplot(3,2,6)
hold on
for i=1:length(recstats)
    needed=find(~isnan(recstats(i).newbaselineSD));
    groupnum=find(strcmp(recstats(i).anaesth,anaesthgroups));
    colorka=colororder{1}(groupnum,:);
    if ~isempty(needed)
        plot(newbaselineSDtimevector(needed),recstats(i).newbaselineSD(needed),'Color',colorka,'LineWidth',2);
    end
end
ylim([0 1.2])
xlabel('time from the start of the recording (min)')
ylabel('relative baseline SD')
saveas(gcf,[dirs.figuresdir,'recording_stats'],'pdf')
saveas(gcf,[dirs.figuresdir,'recording_stats'],'jpg')