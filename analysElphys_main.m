%%
close all
clear all
projectnames={'CB1elphys','InVivo','Persistent-ChRstim','persistent firing','bleb recording'};
% projectnum=3;

% [projectnum,ok] = listdlg('ListString',projectnames,'ListSize',[300 600]); % az XLS file alapján kiválasztjuk, hogy melyik file összes mérésén szeretnénk végigmenni
projectdata.owbridge=0;
projectdata.owbridge=0;
projectdata.owevent=0;
h = aE_projectselector(projectnames);
uiwait(h);
projectnum=projectdata.projectnum;


% alapadatok


if projectnum==1;
    overwrite=0;
    locations=marcicucca_locations;
    dirs.basedir=[locations.tgtardir,'ANALYSISdata/marci/Human_rosehip/CB1elphys/'];
    dirs.rawexporteddir=[dirs.basedir,'Exported_raw/'];
    dirs.bridgeddir=[dirs.basedir,'Bridged_stim/'];
    dirs.eventdir=[dirs.basedir,'Events/'];
    dirs.eventparaleldir=[dirs.basedir,'Events/paralel/'];
    dirs.onlyAPeventdir=[dirs.basedir,'Events_onlyAP/'];
    dirs.grpupedeventdir=[dirs.basedir,'Events_grouped/'];
    dirs.stimepochdir=[dirs.basedir,'Stimepochs/'];
    dirs.figuresdir=[dirs.basedir,'Figures/'];
    amplifier='HEKA';
    xlsdata=aE_readxls([dirs.basedir,'cb1elphys.xls']);
elseif projectnum==2;
    overwrite=0;
    locations=marcicucca_locations;
    dirs.basedir=[locations.tgtardir,'ANALYSISdata/marci/_persistent/_InVivo/'];
    %     dirs.rawdir=[locations.tgtardir,'AXONdata/'];
    dirs.bridgeddir=[dirs.basedir,'Bridged_stim/'];
    dirs.rawexporteddir=[dirs.basedir,'Exported_raw/'];
    dirs.eventdir=[dirs.basedir,'Events/'];
    dirs.eventparaleldir=[dirs.basedir,'Events/paralel/'];
    dirs.statedir=[dirs.basedir,'State/'];
    % dirs.onlyAPeventdir=[dirs.basedir,'Events_onlyAP/'];
    % dirs.grpupedeventdir=[dirs.basedir,'Events_grouped/'];
    % dirs.stimepochdir=[dirs.basedir,'Stimepochs/'];
    dirs.figuresdir=[dirs.basedir,'Figures/'];
    xlsdata=aE_readxls([dirs.basedir,'invivodata.xls']);
    amplifier='HEKA';
elseif projectnum==3;
    overwrite=0;
    locations=marcicucca_locations;
    dirs.basedir=[locations.tgtardir,'ANALYSISdata/marci/_persistent/_ChRstim/'];
    dirs.rawexporteddir=[dirs.basedir,'Exported_raw/'];
    dirs.bridgeddir=[dirs.basedir,'Bridged_stim/'];
    dirs.eventdir=[dirs.basedir,'Events/'];
    dirs.eventparaleldir=[dirs.basedir,'Events/paralel/'];
    %     dirs.onlyAPeventdir=[dirs.basedir,'Events_onlyAP/'];
    %     dirs.grpupedeventdir=[dirs.basedir,'Events_grouped/'];
    %     dirs.stimepochdir=[dirs.basedir,'Stimepochs/'];
    dirs.figuresdir=[dirs.basedir,'Figures/'];
    amplifier='HEKA';
    xlsdata=aE_readxls([dirs.basedir,'ChRstimdata_windows.xls']);
elseif projectnum==4
    overwrite=0;
    locations=marcicucca_locations;
    dirs.basedir=[locations.EMdir,'ANALYSISdata/marci/_persistent/'];
    dirs.rawexporteddir=[dirs.basedir,'Exported_raw/'];
    dirs.bridgeddir=[dirs.basedir,'Bridged_stim/'];
    dirs.eventdir=[dirs.basedir,'Events/'];
    dirs.eventparaleldir=[dirs.basedir,'Events/paralel/'];
    dirs.onlyAPeventdir=[dirs.basedir,'Events_onlyAP/'];
    dirs.grpupedeventdir=[dirs.basedir,'Events_grouped/'];
    dirs.stimepochdir=[dirs.basedir,'Stimepochs/'];
    dirs.figuresdir=[dirs.basedir,'figures/'];
    amplifier='HEKA';
    xlsdata=aE_readxls([dirs.basedir,'persistentdata_windows.xls']);
    dirs.v0distdir=[dirs.basedir,'v0_dist/'];
elseif projectnum==5
    overwrite=0;
    locations=marcicucca_locations;
    dirs.basedir=[locations.tgtardir,'ANALYSISdata/marci/_persistent/_BlebRecording/'];
    dirs.rawexporteddir=[dirs.basedir,'Exported_raw/'];
    dirs.bridgeddir=[dirs.basedir,'Bridged_stim/'];
    dirs.eventdir=[dirs.basedir,'Events/'];
    dirs.eventparaleldir=[dirs.basedir,'Events/paralel/'];
    %     dirs.onlyAPeventdir=[dirs.basedir,'Events_onlyAP/'];
    %     dirs.grpupedeventdir=[dirs.basedir,'Events_grouped/'];
    %     dirs.stimepochdir=[dirs.basedir,'Stimepochs/'];
    dirs.figuresdir=[dirs.basedir,'Figures/'];
    amplifier='HEKA';
    xlsdata=aE_readxls([dirs.basedir,'blebdata_windows.xls']);
end
%%
if strcmp(amplifier,'AXON')
    aE_exportAXONdata(dirs,xlsdata,projectdata.owexport)
else
    %% export Raw data from HEKA
    aE_exportrawHEKAdata(dirs,xlsdata,projectdata.owexport)
    %% generate PGF data, bridge balancing
    overwrite=1;
    aE_generatePGF_bridge_HEKA(dirs,xlsdata,projectdata.owbridge)
end

%% finding events
valtozok.overwrite=projectdata.owevent;
valtozok.plotit=0;
valtozok.threshholdaveragetime=30;%s
% valtozok.mindvpdt=1;
valtozok.minampl=.0001;
valtozok.apampl=.005;
% valtozok.maxfilteredhw=.01;
valtozok.maxdecayriseratio=5;
valtozok.maxrisetime=.01;
valtozok.minaphw=.0001;
valtozok.maxapwidth=.005;
valtozok.filtermovingtime=.0025;
valtozok.diffmovingt=.0005;
valtozok.steptime=.0005; %s
valtozok.eventminsdval=3;
valtozok.apthreshval=10;
parallelcount=4;
aE_findevents(valtozok,dirs,parallelcount,xlsdata)
paralleldata.count=NaN;
while isnan(paralleldata.count) | ~isempty(paralleldata.files)
    paralleldata.files=dir(dirs.eventparaleldir);
    paralleldata.files([paralleldata.files.isdir])=[];
    pause(3)
    paralleldata.prevcount=paralleldata.count;
    paralleldata.count=length(paralleldata.files);
    if paralleldata.prevcount~=paralleldata.count
        disp(['waiting for ',num2str(paralleldata.count),' eventfinding scripts to finish'])
    end
end

%% UP-DOWN transitions
if projectnum==2
    valtozok=struct;
    valtozok.overwrite=0;
    valtozok.segmentlength=5;
    valtozok.plotthestuff=0;
    aE_UP_DOWN_detect(valtozok,dirs,xlsdata);
    
    valtozok=struct;
    valtozok.overwrite=1;
    valtozok.segmentlength=5;
    valtozok.plotthestuff=1;
    aE_UP_DOWN_detect_field(valtozok,dirs,xlsdata);
end

%% plot state transitions
timeborders=[69070, 69510];
timeborders=[0, 76500];
timeborders=[0, inf];
minimumstateduration=0;
[Selection,ok] = listdlg('ListString',{xlsdata.ID},'ListSize',[300 600]);
ID=xlsdata(Selection).ID;
load([dirs.eventdir,ID],'eventdata');
load([dirs.bridgeddir,ID],'bridgeddata','stimdata');
load([dirs.statedir,ID],'statedata');
statedata.UP=statedata.UP([statedata.UP.onsett]>timeborders(1)&[statedata.UP.onsett]<timeborders(2));
statedata.DOWN=statedata.DOWN([statedata.DOWN.onsett]>timeborders(1)&[statedata.DOWN.onsett]<timeborders(2));
% get transition data
timebefore=1;
timeafter=1;
apdata=eventdata(strcmp({eventdata.type},'AP'));
prevsweepnum=0;
for api=1:length(apdata)
    sweepnum=apdata(api).sweepnum;
    if sweepnum~=prevsweepnum
        y=bridgeddata(sweepnum).y;
        si=bridgeddata(sweepnum).si;
        yfiltered=moving(y,5);
        dyfiltered=diff(yfiltered)/si;
        prevsweepnum=sweepnum;
        stepback=round(.0002/si);
    end
    onseth=apdata(api).onseth;
    maxh=apdata(api).maxh;
    threshh=maxh-stepback;
    while dyfiltered(threshh)>5 & threshh>2
        threshh=threshh-1;
    end
    threshv=yfiltered(threshh);
    apdata(api).threshv=threshv;
    apdata(api).APamplitude=apdata(api).maxval-apdata(api).threshv;
    %     if threshv>-.04% & threshv<-.04
    %         figure(1)
    %         clf
    %         hold on
    %         plot(yfiltered,'r-')
    %         plot(y,'k-')
    %
    %         plot(maxh,yfiltered(maxh),'ro')
    %         plot(threshh,yfiltered(threshh),'ko')
    %         plot(onseth,yfiltered(onseth),'kx')
    %         pause
    %     end
end
figure(1)
clf
% hist([apdata.threshv],100)
plot([apdata.threshv],[apdata.APamplitude],'ko')
[x,~]=ginput(1);
aapdata=apdata([apdata.threshv]<x);
sapdata=apdata([apdata.threshv]>x);
epdata=eventdata(strcmp({eventdata.type},'ep'));
ipdata=eventdata(strcmp({eventdata.type},'ip'));
statestocheck={'UP','DOWN'};
%
for statei=1:length(statestocheck)
    statedatanow=statedata.(statestocheck{statei});
    statedatanow=statedatanow([statedatanow.duration]>minimumstateduration);
    transitiondata=struct;
    for i=1:length(statedatanow)
        sweepnum=statedatanow(i).sweepnum;
        transitiont=statedatanow(i).onsett;
        transitionh=statedatanow(i).onseth;
        si=bridgeddata(sweepnum).si;
        stepback=round(timebefore/si);
        stepforward=round(timeafter/si);
        y=bridgeddata(sweepnum).y;
        if transitionh>stepback & length(y)>transitionh+stepforward
            if isempty(fieldnames(transitiondata))
                NEXT=1;%% recording statistics
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
                figure(1)
                clf
                subplot(3,1,1)
                [nall,xbin]=hist([recstats.recordinglength]/60,[2.5:5:62.5]);
                [nketamine,~]=hist([recstats(strcmp({recstats.anaesth},'ketamine xylazine')).recordinglength]/60,[2.5:5:62.5]);
                [nchloral,~]=hist([recstats(strcmp({recstats.anaesth},'chloral hydrate')).recordinglength]/60,[2.5:5:62.5]);
                [hAxes,hBar,hLine] = plotyy(xbin,[nketamine;nchloral]',xbin,cumsum([nketamine;nchloral]'),'bar','plot');
                set(hBar,'BarLayout','stacked')
                colororder=get(hAxes,'colororder');
                set(hBar(2),'FaceColor',colororder{1}(2,:));
                set(hLine,'LineWidth',2)
                legend({'ketamine xylazine','chloral hydrate'})
                set(gca,'Xtick',[0:5:60])
                xlabel('min')
                ylabel('# of cells')
                title('recording length')
                subplot(3,1,2)
                [nall,xbin]=hist([recstats.RS]/10^6,[5:10:95]);
                [nketamine,~]=hist([recstats(strcmp({recstats.anaesth},'ketamine xylazine')).RS]/10^6,[5:10:95]);
                [nchloral,~]=hist([recstats(strcmp({recstats.anaesth},'chloral hydrate')).RS]/10^6,[5:10:95]);
                [hAxes,hBar,hLine] = plotyy(xbin,[nketamine;nchloral]',xbin,cumsum([nketamine;nchloral]'),'bar','plot');
                set(hBar,'BarLayout','stacked')
                colororder=get(hAxes,'colororder');
                set(hBar(2),'FaceColor',colororder{1}(2,:));
                set(hLine,'LineWidth',2)
                legend({'ketamine xylazine','chloral hydrate'})
                set(gca,'Xtick',[0:10:90])
                xlabel('MOhm')
                ylabel('# of cells')
                title('median RS')
                subplot(3,1,3)
                hold on
                for i=1:length(recstats)
                    needed=find(~isnan(recstats(i).newbaselineSD));
                    if strcmp(recstats(i).anaesth,'ketamine xylazine')
                        colorka=colororder{1}(1,:);
                    elseif strcmp(recstats(i).anaesth,'chloral hydrate')
                        colorka=colororder{1}(2,:);
                    end
                    if ~isempty(needed)
                        plot(newbaselineSDtimevector(needed),recstats(i).newbaselineSD(needed),'Color',colorka,'LineWidth',2);
                    end
                end
                ylim([0 1.2])
                xlabel('time (min)')
                ylabel('relative baseline SD')
            else
                NEXT=length(transitiondata)+1;
            end
            neededAPs=find([apdata.maxtime]>transitiont-timebefore &[apdata.maxtime]<transitiont+timeafter);
            neededaAPs=find([aapdata.maxtime]>transitiont-timebefore &[aapdata.maxtime]<transitiont+timeafter);
            neededsAPs=find([sapdata.maxtime]>transitiont-timebefore &[sapdata.maxtime]<transitiont+timeafter);
            neededeps=find([epdata.maxtime]>transitiont-timebefore &[epdata.maxtime]<transitiont+timeafter);
            neededips=find([ipdata.maxtime]>transitiont-timebefore &[ipdata.maxtime]<transitiont+timeafter);
            transitiondata(NEXT).transitiont=transitiont;
            transitiondata(NEXT).time=[-stepback:stepforward]'*si;
            transitiondata(NEXT).y=y(transitionh-stepback:transitionh+stepforward)';
            
            if ~isempty(neededAPs)
                transitiondata(NEXT).APtimes=[apdata(neededAPs).maxtime]-transitiont;
                transitiondata(NEXT).APnum=length(transitiondata(NEXT).APtimes);
            else
                transitiondata(NEXT).APtimes=[];
                transitiondata(NEXT).APnum=0;
            end
            if ~isempty(neededaAPs)
                transitiondata(NEXT).aAPtimes=[aapdata(neededaAPs).maxtime]-transitiont;
                transitiondata(NEXT).aAPnum=length(transitiondata(NEXT).aAPtimes);
            else
                transitiondata(NEXT).aAPtimes=[];
                transitiondata(NEXT).aAPnum=0;
            end
            if ~isempty(neededsAPs)
                transitiondata(NEXT).sAPtimes=[sapdata(neededsAPs).maxtime]-transitiont;
                transitiondata(NEXT).sAPnum=length(transitiondata(NEXT).sAPtimes);
            else
                transitiondata(NEXT).sAPtimes=[];
                transitiondata(NEXT).sAPnum=0;
            end
            if ~isempty(neededeps)
                transitiondata(NEXT).eptimes=[epdata(neededeps).maxtime]-transitiont;
                transitiondata(NEXT).epnum=length(transitiondata(NEXT).eptimes);
            else
                transitiondata(NEXT).eptimes=[];
                transitiondata(NEXT).epnum=0;
            end
            if ~isempty(neededips)
                transitiondata(NEXT).iptimes=[ipdata(neededips).maxtime]-transitiont;
                transitiondata(NEXT).ipnum=length(transitiondata(NEXT).iptimes);
            else
                transitiondata(NEXT).iptimes=[];
                transitiondata(NEXT).ipnum=0;
            end
        end
    end
    Transitiondata.(statestocheck{statei})=transitiondata;
end
% plot
figure(1)
clf
subplot(4,2,1)
needed=find([Transitiondata.UP.aAPnum]>0);
plot([Transitiondata.UP(needed).time],[Transitiondata.UP(needed).y]);
xlim([-timebefore,timeafter])
title('DOWN to UP transition')
subplot(4,2,2)
needed=find([Transitiondata.DOWN.aAPnum]>0);
plot([Transitiondata.DOWN(needed).time],[Transitiondata.DOWN(needed).y]);
xlim([-timebefore,timeafter])
title('UP to DOWN transition')
subplot(4,2,3)
hist([Transitiondata.UP.sAPtimes],[-timebefore:0.05:timeafter])
xlim([-timebefore,timeafter])
ylabel('somatic APnum');
subplot(4,2,4)
hist([Transitiondata.DOWN.sAPtimes],[-timebefore:0.05:timeafter])
xlim([-timebefore,timeafter])
ylabel('somatic AP num');
subplot(4,2,5)
hist([Transitiondata.UP.aAPtimes],[-timebefore:0.05:timeafter])
xlim([-timebefore,timeafter])
ylabel('axonal APnum');
subplot(4,2,6)
hist([Transitiondata.DOWN.aAPtimes],[-timebefore:0.05:timeafter])
xlim([-timebefore,timeafter])
ylabel('axonal AP num');
subplot(4,2,7)
hist([Transitiondata.UP.eptimes],[-timebefore:0.05:timeafter])
xlim([-timebefore,timeafter])
ylabel('EPnum');
subplot(4,2,8)
hist([Transitiondata.DOWN.eptimes],[-timebefore:0.05:timeafter])
xlim([-timebefore,timeafter])
ylabel('EPnum');

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
figure(1)
clf
subplot(4,1,1)
[nall,xbinn]=hist([recstats.recordinglength]/60,[2.5:5:62.5]);
[RSall,xbinRS]=hist([recstats.RS]/10^6,[5:10:95]);
xbinZ=[5:10:155];
histstruct=struct;
for groupi=1:length(anaesthgroups)
    [histstruct(groupi).n,~]=hist([recstats(strcmp({recstats.anaesth},anaesthgroups{groupi})).recordinglength]/60,xbinn);
    histstruct(groupi).n=histstruct(groupi).n';
    [histstruct(groupi).RS,~]=hist([recstats(strcmp({recstats.anaesth},anaesthgroups{groupi})).RS]/10^6,xbinRS);
    histstruct(groupi).RS=histstruct(groupi).RS';
    histstruct(groupi).Z=[];
    for xlsi=1:length(xlsdata)
        if isnumeric( xlsdata(xlsi).locationz) & strcmp(xlsdata(xlsi).anaesthesia,anaesthgroups{groupi})
            histstruct(groupi).Z=[histstruct(groupi).Z, xlsdata(xlsi).locationz];
        end
    end
    [histstruct(groupi).Zhist,~]=hist(histstruct(groupi).Z,xbinZ);
    histstruct(groupi).Zhist=histstruct(groupi).Zhist';
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


subplot(4,1,2)
[hAxes,hBar,hLine] = plotyy(xbinRS,[histstruct.RS],xbinRS,cumsum([histstruct.RS]),'bar','plot');
set(hBar,'BarLayout','stacked')
set(gca,'box','off')
% ylim([0,max(sum([histstruct.RS]'))]);
set(hAxes(1), 'YLimMode', 'Auto');
set(hAxes(1),'Ytickmode','auto');
set(hAxes(2),'Ytickmode','auto');
set(hAxes(2), 'YLimMode', 'Auto');
set(hAxes(2),'XTick',[],'XTickLabel',[]);
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


subplot(4,1,3)
[hAxes,hBar,hLine] = plotyy(xbinZ,[histstruct.Zhist],xbinZ,cumsum([histstruct.Zhist]),'bar','plot');
set(hBar,'BarLayout','stacked')
set(gca,'box','off')
% ylim([0,max(sum([histstruct.RS]'))]);
set(hAxes(1), 'YLimMode', 'Auto');
set(hAxes(1),'Ytickmode','auto');
set(hAxes(2),'Ytickmode','auto');
set(hAxes(2), 'YLimMode', 'Auto');
set(hAxes(2),'XTick',[],'XTickLabel',[]);
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


subplot(4,1,4)
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
%%
if projectnum==4
    %% defining stimepochs and spike clusters
    valtozok_stimepochs.overwrite=projectdata.owstimepoch;
    valtozok_stimepochs.histbins=[0:.001:.5];
    valtozok_stimepochs.sdtimesval_persistentgroup=4;
    valtozok_stimepochs.maxskewness_persistentgroup=2;
    valtozok_stimepochs.maxisiwithingroup=.1;
    valtozok_stimepochs.minapnum_persistentgroup=1;
    valtozok_stimepochs.minAPtimediff_stimepoch=10;
    valtozok_stimepochs.minAPnum_stimepoch=20;
    valtozok_stimepochs.plotpersistentgroups=0;
    persistent_definestimepoch(dirs,valtozok_stimepochs)
    
    %% determining TAU and V0
    valtozok_tau.timeconstanthossz=.100; %s
    valtozok_tau.stepstocheck=20;%*2
    valtozok_states.overwrite=0;
    valtozok_states.ordfiltorder=.2;
    valtozok_states.ordfiltlength=.02; %s
    files=dir(dirs.bridgeddir);
    files([files.isdir])=[];
    for fnum=1:length(files)
        progressbar(fnum/length(files),[],[])
        fname=files(fnum).name;
        a=dir([dirs.v0distdir,fname]);
        if isempty(a) | valtozok_states.overwrite==1
            taudata=persistent_gettau(valtozok_tau,dirs,fname);
            persistent_getv0dist(dirs,valtozok_states,valtozok_tau,taudata,fname);
        end
    end
end
if  projectnum==5
    %% checking somatically triggered axonal events
    timebefore=.003;
    timeafter=.005;
    cutofffreq=[3000];
    filterorder=1;
    for i=1:length(xlsdata)
        if xlsdata(i).startT>xlsdata(i).endT
            xlsdata(i).endTmodified=xlsdata(i).endT+86400;
        else
            xlsdata(i).endTmodified=xlsdata(i).endT;
        end
    end
    blebidxes=find([xlsdata.axonal]);
    fieldidxes=find([xlsdata.field]);
    icdxes=find(~[xlsdata.field]&~[xlsdata.axonal]);
    [icSelection,ok] = listdlg('ListString',{xlsdata(icdxes).ID},'ListSize',[300 600]); % az XLS file alapján kiválasztjuk, hogy melyik file összes mérésén szeretnénk végigmenni
    for ici=1:length(icSelection)
        icxlsidx=icdxes(icSelection(ici));
        blebSelection=zeros(size(blebidxes));
        for blebi=1:length(blebidxes)
            blebidx=blebidxes(blebi);
            if strcmp(xlsdata(icxlsidx).HEKAfname,xlsdata(blebidx).HEKAfname) & xlsdata(blebidx).startT<xlsdata(icxlsidx).endTmodified & xlsdata(icxlsidx).startT<xlsdata(blebidx).endTmodified
                blebSelection(blebi)=1;
            end
        end
        blebSelection=blebidxes(find(blebSelection));
        for blebi=1:length(blebSelection)
            blebxlsidx=blebSelection(blebi);
            disp(xlsdata(blebxlsidx).ID)
            ic=struct;
            bleb=struct;
            ic=load([dirs.bridgeddir,xlsdata(icxlsidx).ID]);
            load([dirs.eventdir,xlsdata(icxlsidx).ID],'eventdata');
            ic.eventdata=eventdata;
            bleb=load([dirs.bridgeddir,xlsdata(blebxlsidx).ID]);
            for blebsweepnum=1:length(bleb.bridgeddata) %filtering
                if strcmp(xlsdata(blebxlsidx).ID,'1702083rm_2_20_1')
                    bleb.bridgeddata(blebsweepnum).realtime=bleb.bridgeddata(blebsweepnum).realtime+24*3600;
                end
                si=bleb.bridgeddata(blebsweepnum).si;
                ninq=.5/si;
                if length(cutofffreq)==1
                    [b,a] = butter(filterorder,cutofffreq/ninq,'low');
                else
                    [b,a] = butter(filterorder,cutofffreq/ninq,'bandpass');
                end
                bleb.bridgeddata(blebsweepnum).y_filt=filter(b,a,bleb.bridgeddata(blebsweepnum).y);
            end
            apidxes=find(strcmp({ic.eventdata.type},'AP'));
            APdata=struct;
            NEXT=0;
            for api=1:length(apidxes)
                eventidx=apidxes(api);
                si=round(ic.eventdata(eventidx).si*10^6)/10^6;
                icsweepnum=ic.eventdata(eventidx).sweepnum;
                maxh=ic.eventdata(eventidx).maxh;
                maxtime=ic.eventdata(eventidx).maxtime;
                stimulated=ic.eventdata(eventidx).stimulated;
                baselineval=ic.eventdata(eventidx).baselineval;
                amplitude=ic.eventdata(eventidx).amplitude;
                halfwidth=ic.eventdata(eventidx).halfwidth;
                blebsweepnum=find([bleb.bridgeddata.realtime]==ic.bridgeddata(icsweepnum).realtime);
                stepback=round(timebefore/si);
                stepforward=round(timeafter/si);
                
                if ~isempty(blebsweepnum) && maxh>stepback && stepforward+maxh<length(ic.bridgeddata(icsweepnum).y)
                    NEXT=NEXT+1;
                    APdata(NEXT).maxtime=maxtime;
                    APdata(NEXT).si=si;
                    APdata(NEXT).time=[-stepback:stepforward]'*si*1000;
                    APdata(NEXT).ic_y=ic.bridgeddata(icsweepnum).y(maxh-stepback:maxh+stepforward)';
                    APdata(NEXT).bleb_y=bleb.bridgeddata(blebsweepnum).y(maxh-stepback:maxh+stepforward)';
                    APdata(NEXT).bleb_y_filt=bleb.bridgeddata(blebsweepnum).y_filt(maxh-stepback:maxh+stepforward)';
                    APdata(NEXT).bleb_y_baseline=mean(APdata(NEXT).bleb_y(1:stepback));
                    APdata(NEXT).bleb_channellabel=bleb.lightdata(blebsweepnum).channellabel;
                    APdata(NEXT).bleb_preamplnum=bleb.lightdata(blebsweepnum).preamplnum;
                    APdata(NEXT).bleb_Amplifiermode=bleb.lightdata(blebsweepnum).Amplifiermode;
                    APdata(NEXT).stimulated=stimulated;
                    APdata(NEXT).baselineval=median(APdata(NEXT).ic_y(1:stepback));
                    APdata(NEXT).amplitude=amplitude;
                    APdata(NEXT).halfwidth=halfwidth;
                end
            end
            si=max([APdata.si]);
            for NEXT=1:length(APdata) %resampling
                if APdata(NEXT).si<si
                    APdata(NEXT).time=resample(APdata(NEXT).time,APdata(NEXT).si*10^6,si*10^6);
                    APdata(NEXT).ic_y=resample(APdata(NEXT).ic_y,APdata(NEXT).si*10^6,si*10^6);
                    APdata(NEXT).bleb_y=resample(APdata(NEXT).bleb_y,APdata(NEXT).si*10^6,si*10^6);
                    APdata(NEXT).bleb_y_filt=resample(APdata(NEXT).bleb_y_filt,APdata(NEXT).si*10^6,si*10^6);
                    APdata(NEXT).si=si;
                end
            end
        end
    end
    APdataOriginal=APdata;
    %% visualize bleb data
    APdata=APdataOriginal;
    recordingmode='V-Clamp';
    baselinepercentilerange=[.01 .99];
    baselinestdmin=3;
    % detect amplitudes
    %      figure(1)
    %     subplot(4,4,i*4+2);
    %     [x,~]=ginput(2);
    x=[-2,2];
    baselinex=[-2, 0];
    for api=1:length(APdata)
        bleb_y_baseline=APdata(api).bleb_y_baseline;
        time=APdata(api).time;
        bleb_y=APdata(api).bleb_y;
        bleb_y_filt=APdata(api).bleb_y_filt;
        ic_y=APdata(api).ic_y;
        idxes=find(time>min(x)&time<max(x));
        baselineidxes=find(time>min(baselinex)&time<max(baselinex));
        if strcmp(recordingmode,'V-Clamp')
            [peakv,peakh]=max(bleb_y_filt(idxes));
            peakv=(peakv-bleb_y_baseline);
            peakh=peakh+idxes(1);
            [baselinepeakv,~]=min(bleb_y_filt(baselineidxes));
            baselinepeakv=-(baselinepeakv-bleb_y_baseline);
        elseif strcmp(recordingmode,'C-Clamp')
            [peakv,peakh]=min(bleb_y_filt(idxes));
            peakv=-(peakv-bleb_y_baseline);
            peakh=peakh+idxes(1);
            [baselinepeakv,~]=max(bleb_y_filt(baselineidxes));
            baselinepeakv=(baselinepeakv-bleb_y_baseline);
        end
        APdata(api).bleb_amplitude=peakv;
        APdata(api).bleb_amplitude_t=time(peakh);
        APdata(api).bleb_baselineamplitude=baselinepeakv;
    end
    
    baselineamplitudes=sort([APdata(strcmp({APdata.bleb_Amplifiermode},recordingmode)).bleb_baselineamplitude]);
    
    [nbase,xbincenters]=hist(baselineamplitudes,100);
    baselineamplitudes=baselineamplitudes(round(baselinepercentilerange(1)*length(baselineamplitudes)):round(baselinepercentilerange(2)*length(baselineamplitudes)));
    [nbase2,~]=hist(baselineamplitudes,xbincenters);
    minamplitude=mean(baselineamplitudes)+std(baselineamplitudes)*baselinestdmin;
    
    todel=find([APdata.bleb_amplitude]<minamplitude&strcmp({APdata.bleb_Amplifiermode},recordingmode) );
    
    APdata(todel)=[];
    
    close all
    figure(33)
    clf
    bar(xbincenters,nbase,'r')
    hold on
    bar(xbincenters,nbase2,'b')
    plot([minamplitude,minamplitude],[0, max(nbase)],'r-','LineWidth',3)
    figure(1)
    clf
    hold on
    figure(2)
    clf
    hold on
    for i =0:3
        figure(1)
        needed=(strcmp({APdata.bleb_Amplifiermode},recordingmode));
        if i==0 % persistent spikelet
            cim='axonal spikelet';
            needed=find(needed&[APdata.stimulated]==0 & [APdata.amplitude]<.06 & [APdata.halfwidth]>300*10^-6); %
        elseif i==1 % persistent spike
            cim='axonal spike';
            needed=find(needed&[APdata.stimulated]==0 & [APdata.amplitude]>.07 &  [APdata.baselineval]<-.00); %& [APdata.halfwidth]>300*10^-6
        elseif i==2 % somatic spike
            cim='somatic spike - hyperpolarized v0';
            needed=find(needed&[APdata.stimulated]==1 & [APdata.amplitude]>.07&  [APdata.baselineval]<-.07);
        elseif i==3 % somatic spike
            cim='somatic spike - depolarized v0';
            needed=find(needed&[APdata.stimulated]==1 & [APdata.amplitude]>.07&  [APdata.baselineval]>-.05);
        end
        figure(1)
        hsom(i+1)=subplot(4,4,i*4+1);
        hold on
        plot([APdata(needed).time],[APdata(needed).ic_y])
        plot(mean([APdata(needed).time],2),mean([APdata(needed).ic_y],2),'k-','LineWidth',3)
        title(cim)
        hax(i+1)=subplot(4,4,i*4+2);
        hold on
        plot([APdata(needed).time],bsxfun(@(a,b)a-b,[APdata(needed).bleb_y_filt],[APdata(needed).bleb_y_baseline]))
        plot(mean([APdata(needed).time],2),mean(bsxfun(@(a,b)a-b,[APdata(needed).bleb_y_filt],[APdata(needed).bleb_y_baseline]),2),'k-','LineWidth',3)
        if isfield(APdata,'bleb_amplitude')
            [~,amplbins]=hist([APdata(strcmp({APdata.bleb_Amplifiermode},recordingmode)).bleb_amplitude],100);
            [~,ampltbins]=hist([APdata(strcmp({APdata.bleb_Amplifiermode},recordingmode)).bleb_amplitude_t],100);
            subplot(4,4,i*4+3);
            hist([APdata(needed).bleb_amplitude],amplbins)
            xlabel('peak amplitude')
            subplot(4,4,i*4+4);
            hist([APdata(needed).bleb_amplitude_t],ampltbins)
            xlabel('peak latency')
        end
        linkaxes(hsom,'xy')
        linkaxes(hax,'xy')
        subplot(4,4,i*4+2);
        axis tight
        subplot(4,4,i*4+1);
        axis tight
        if isfield(APdata,'bleb_amplitude')
            figure(2)
            if i<2
                colorka='rx';
            else
                colorka='ko';
            end
            subplot(3,1,1)
            hold on
            plot([APdata(needed).maxtime],[APdata(needed).bleb_amplitude],colorka)
            ylabel('peak amplitude')
            subplot(3,1,2)
            hold on
            plot([APdata(needed).maxtime],[APdata(needed).bleb_amplitude_t],colorka)
            ylabel('peak latency')
            xlabel('time')
            subplot(3,1,3)
            hold on
            plot([APdata(needed).bleb_amplitude_t],[APdata(needed).bleb_amplitude],colorka)
            xlabel('peak latency')
            ylabel('peak amplitude')
        end
    end
    
end
% return
% %% uj cucc
%
% files=dir(dirs.bridgeddir);
% files([files.isdir])=[];
% sweepdata=struct;
% progressbar('preloading files')
% for fi=1:length(files)
%     load([dirs.bridgeddir,files(fi).name],'lightdata');
%     if isempty(fieldnames(sweepdata))
%         sweepdata=lightdata;
%     else
%         sweepdata=[sweepdata,lightdata];
%     end
%     progressbar(fi/length(files))
% end
%
%
% %% puffnlightstim
% aE_persistent_puffnlightstim %this script plots puffing and light stim experiments.. ap waveforms and onsets are analysed - base values are needed from this main script
%% check electrotonic and chemical connectivity
valtozok.plot.dpi=150;
valtozok.plot.xcm=20;
valtozok.plot.ycm=14;
valtozok.plot.betumeret=14;
valtozok.plot.betutipus='Arial';
valtozok.plot.axesvastagsag=2;
valtozok.plot.xinch=valtozok.plot.xcm/2.54;
valtozok.plot.yinch=valtozok.plot.ycm/2.54;
valtozok.plot.xsize=valtozok.plot.dpi*valtozok.plot.xinch;
valtozok.plot.ysize=valtozok.plot.dpi*valtozok.plot.yinch;


valtozok.gj_baselinelength=.010;
valtozok.gj_baselinelengthend=.08;
valtozok.gj_minlinelength=.05;
valtozok.gj_mincurrampl=-10*10^-12;

valtozok.noAPbeforetheevent=1; %s
valtozok.noAPaftertheevent=.05; %s
valtozok.pairedpulseneeded=0; %boolean
valtozok.pairedpulsedelay=NaN;%.06; %s
valtozok.pairedpulsejitter=.05; %s
valtozok.baselinelength=0.025; %s
valtozok.psplength=.15; %s
valtozok.filterorder=3;
valtozok.cutofffreq=1000;
valtozok.drugwashintime=120;
valtozok.maxy0baselinedifference=.0005;
valtozok.discardpostsweepswithap=1;
valtozok.postrecordingmode='C-Clamp';%'C-Clamp' or 'V-Clamp' or 'any'
valtozok.prerecordingmode='C-Clamp';%'C-Clamp' or 'V-Clamp' or 'any'

aE_checkGJandChemicalSynapse(valtozok,xlsdata,dirs)
%% plotting IV

valtozok.plot.betumeret=8;
valtozok.plot.axesvastagsag=2;
valtozok.plot.xcm=20;
valtozok.plot.ycm=14;
valtozok.plot.dpi=150;
xinch=valtozok.plot.xcm/2.54;
yinch=valtozok.plot.ycm/2.54;
valtozok.plot.xsize=valtozok.plot.dpi*xinch;
valtozok.plot.ysize=valtozok.plot.dpi*yinch;


dothesecond=zeros(size(xlsdata));

for prenum=1:length(xlsdata)
    %     pause
    fname=[xlsdata(prenum).ID,'.mat'];
    HEKAfname=xlsdata(prenum).HEKAfname;
    if any(strfind(HEKAfname,','))
        comma=strfind(HEKAfname,',');
        HEKAfname=HEKAfname(1:comma(1)-1);
    end
    load([locations.tgtardir,'MATLABdata/IV/',xlsdata(prenum).setup,'/',HEKAfname,'.mat']);
    
    gsc=xlsdata(prenum).G_S_C;
    hyps=strfind(gsc,'_');
    commas=strfind(gsc,',');
    g=num2str(gsc(1:hyps(1)-1));

    commas=[hyps(1),commas,hyps(2)];
    for is=1:length(commas)-1
        s=num2str(gsc(commas(is)+1:commas(is+1)-1));
        c=num2str(gsc(hyps(2)+1:end));
        iv=iv.(['g',g,'_s',s,'_c',c]);
        si=mode(diff(iv.time));
        [b,a]=butter(3,15000/(1/mode(diff(iv.time)))/2,'low');
        [bb,aa]=butter(3,1000/(1/mode(diff(iv.time)))/2,'low');
            
        

        figure(33)
        clf
        hold on
        plot(iv.time,filtfilt(bb,aa,iv.v1),'k-','LineWidth',2)
        plot(iv.time,filtfilt(b,a,iv.(['v',num2str(iv.sweepnum-(sweeplevonas))])),'k-','LineWidth',2);
        axis tight
        ylim([-.13 .050])
        xlim([0 1])
        %     ylimitek(:,i)=get(gca,'Ylim');
        set(gca,'LineWidth',valtozok.plot.axesvastagsag,'FontSize',valtozok.plot.betumeret,'Position',[1/valtozok.plot.xcm 1/valtozok.plot.ycm 1-2/valtozok.plot.xcm 1-2/valtozok.plot.ycm])
        axis off
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 valtozok.plot.xsize/valtozok.plot.dpi valtozok.plot.ysize/valtozok.plot.dpi])
        print(gcf,[dirs.figuresdir,'/IVs/IV_',xlsdata(prenum).ID,'.jpg'],'-djpeg',['-r',num2str(valtozok.plot.dpi)])
        saveas(gcf,[dirs.figuresdir,'/IVs/IV_',xlsdata(prenum).ID,'.fig'])
    end
end
figure(1)
clf
hold on

plot([.5 .5],[-.050 -.03],'k-','LineWidth',5)
plot([.5 .6],[-.050 -.050],'k-','LineWidth',5)
ylim([-.13 .050])
xlim([0 1])
axis off

set(gca,'Position',[0 0 1 1])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 valtozok.plot.xsize/valtozok.plot.dpi valtozok.plot.ysize/valtozok.plot.dpi])
print(gcf,[dirs.figuresdir,'/IVs/IVscalebar_100ms_20mV.jpg'],'-djpeg',['-r',num2str(valtozok.plot.dpi)])
saveas(gcf,[dirs.figuresdir,'/IVs/IVscalebar_100ms_20mV.fig'])
