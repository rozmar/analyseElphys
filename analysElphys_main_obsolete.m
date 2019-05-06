% everything is moved here from anaysElphys_main.m which are not used
% currently



%%

valtozok=struct;
valtozok.movingvindowsize=3;
valtozok.movingvindowstep=.5; %seconds for downsampling and median filtering
valtozok.timeborders=[0 inf];
valtozok.frequencyrange=[1 4];
valtozok.PSDonfield=true;
valtozok.minsweeplength=valtozok.movingvindowsize/2;
xlsnum=find(strcmp({xlsdata.ID},'1705301rm_1_2_1'));

aE_V0_vs_PSD(dirs,xlsdata,xlsnum,valtozok)

%
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


%%
apersratio=[APdata.aAPnum]./[APdata.sAPnum];
[~,idx]=sort(apersratio,'descend');
APdata=APdata(idx);
% for filei=1:length(files)
%     fname=files(filei).name(1:end-4);
%     load([sortedeventdir,fname],'eventdata');
%     xlsidx=find(strcmp({xlsdata.ID},fname));
%     
%     for xlsfieldi=1:length(fieldek)
%         APdata(filei).(fieldek{xlsfieldi})=xlsdata(xlsidx).(fieldek{xlsfieldi});
%     end
% %     a=dir()
%     APdata(filei).aAPnum=sum([eventdata.axonalAP]);
%     APdata(filei).sAPnum=sum([eventdata.somaticAP]);
% end

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
    cutofffreq=[10000];
    filterorder=1;
    filtertype='gauss';%'butter'
    gaussigma=40;%microsec
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
            disp(['bleb recorded: ', xlsdata(blebxlsidx).ID])
            ic=struct;
            bleb=struct;
            ic=load([dirs.bridgeddir,xlsdata(icxlsidx).ID]);
            load([dirs.eventdir,'sorted/',xlsdata(icxlsidx).ID],'eventdata');
            ic.eventdata=eventdata;
            bleb=load([dirs.bridgeddir,xlsdata(blebxlsidx).ID]);
            for blebsweepnum=1:length(bleb.bridgeddata) %filtering
                if strcmp(xlsdata(blebxlsidx).ID,'1702083rm_2_20_1')
                    bleb.bridgeddata(blebsweepnum).realtime=bleb.bridgeddata(blebsweepnum).realtime+24*3600;
                end
                si=bleb.bridgeddata(blebsweepnum).si;
                if strcmp(filtertype,'butter')
                    ninq=.5/si;
                    if length(cutofffreq)==1
                        [b,a] = butter(filterorder,cutofffreq/ninq,'low');
                    else
                        [b,a] = butter(filterorder,cutofffreq/ninq,'bandpass');
                    end
                    bleb.bridgeddata(blebsweepnum).y_filt=filter(b,a,bleb.bridgeddata(blebsweepnum).y)';
                else
                    movingn=ceil(gaussigma/si*10/1e6);
                    hgauss=fspecial('gaussian',[movingn,1],gaussigma/si/1e6);
                    bleb.bridgeddata(blebsweepnum).y_filt=imfilter(bleb.bridgeddata(blebsweepnum).y',hgauss,'replicate')';
                end
            end
            apidxes=find(strcmp({ic.eventdata.type},'AP'));
            APdata=struct;
            NEXT=0;
            for api=1:length(apidxes)
                eventidx=apidxes(api);
                si=round(ic.eventdata(eventidx).si*10^6)/10^6;
                icsweepnum=ic.eventdata(eventidx).sweepnum;
                maxh=ic.eventdata(eventidx).threshh;%threshh;%maxh;
                maxtime=ic.eventdata(eventidx).maxtime;
                stimulated=ic.eventdata(eventidx).stimulated;
                baselineval=ic.eventdata(eventidx).baselineval;
                amplitude=ic.eventdata(eventidx).APamplitude;
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
                    APdata(NEXT).threshv=eventdata(eventidx).threshv;
                    APdata(NEXT).axonalAP=eventdata(eventidx).axonalAP;
                    APdata(NEXT).somaticAP=eventdata(eventidx).somaticAP;
%                     figure(3)
%                     clf
%                     subplot(3,1,1)
%                     plot(APdata(NEXT).time,APdata(NEXT).ic_y)
%                     subplot(3,1,2)
%                     plot(APdata(NEXT).time,APdata(NEXT).bleb_y)
%                     subplot(3,1,3)
%                     plot(APdata(NEXT).time,APdata(NEXT).bleb_y_filt)
%                     pause
%                     disp('lol')
                end
            end
            si=max([APdata.si]);
            for NEXT=1:length(APdata) %resampling
                if APdata(NEXT).si<si
%                     APdata(NEXT).time=resample(APdata(NEXT).time,APdata(NEXT).si*10^6,si*10^6);
%                     APdata(NEXT).ic_y=resample(APdata(NEXT).ic_y,APdata(NEXT).si*10^6,si*10^6);
%                     APdata(NEXT).bleb_y=resample(APdata(NEXT).bleb_y,APdata(NEXT).si*10^6,si*10^6);
%                     APdata(NEXT).bleb_y_filt=resample(APdata(NEXT).bleb_y_filt,APdata(NEXT).si*10^6,si*10^6);
                    
                    APdata(NEXT).time=downsample(APdata(NEXT).time,round(si/APdata(NEXT).si));
                    APdata(NEXT).ic_y=downsample(APdata(NEXT).ic_y,round(si/APdata(NEXT).si));
                    APdata(NEXT).bleb_y=downsample(APdata(NEXT).bleb_y,round(si/APdata(NEXT).si));
                    APdata(NEXT).bleb_y_filt=downsample(APdata(NEXT).bleb_y_filt,round(si/APdata(NEXT).si));
                    
                    APdata(NEXT).si=si;

                end
            end
        end
    end
    APdataOriginal=APdata;
    %% visualize bleb data
    APdata=APdataOriginal;
    APdata=APdata([APdata.maxtime]>87120 & [APdata.maxtime]<87250 ); % 1702083rm_2_20_1
%     APdata=APdata([APdata.maxtime]>87300 & [APdata.maxtime]<872500 ); % 1702083rm_2_20_1
%      APdata=APdata([APdata.maxtime]>79250 & [APdata.maxtime]<79300 );      %1712103rm
%      APdata=APdata([APdata.maxtime]>79300 & [APdata.maxtime]<79355 ); % 1712103rm
%     APdata=APdata([APdata.maxtime]>79392 & [APdata.maxtime]<79410 ); % 1712103rm
    recordingmode='V-Clamp';
    baselinepercentilerange=[.01 .99];
    baselinestdmin=-3;
    baselinestdmin_for_noise=100;
    minhalfwidth=.1/1000;
    % detect amplitudes
    %      figure(1)
    %     subplot(4,4,i*4+2);
    %     [x,~]=ginput(2);
    x=[0,2];
    baselinex=[-2, 0];
    for api=1:length(APdata)
        bleb_y_baseline=APdata(api).bleb_y_baseline;
        time=APdata(api).time;
        bleb_y=APdata(api).bleb_y;
        bleb_y_filt=APdata(api).bleb_y_filt;
        ic_y=APdata(api).ic_y;
        idxes=find(time>min(x)&time<max(x));
        idxes=find(time>min(x)&time<max(x));
        baselineidxes=find(time>min(baselinex)&time<max(baselinex));
        if strcmp(recordingmode,'V-Clamp')
%             [peakv,peakh]=max(bleb_y_filt(idxes));
%             peakv=(peakv-bleb_y_baseline);
%             peakh=peakh+idxes(1);
%             [baselinepeakv,~]=min(bleb_y_filt(baselineidxes));
%             baselinepeakv=-(baselinepeakv-bleb_y_baseline);
% % % % % % % % %             [peakv,peakh]=min(bleb_y_filt(idxes));
% % % % % % % % %             peakv=(peakv-bleb_y_baseline);
% % % % % % % % %             peakh=peakh+idxes(1);
% % % % % % % % %             [baselinepeakv,~]=max(bleb_y_filt(baselineidxes));
% % % % % % % % %             baselinepeakv=-(baselinepeakv-bleb_y_baseline);
% %             figure(33)
% %             clf
% %             subplot(2,1,1)
% %             plot(APdata(api).ic_y)
% %             subplot(2,1,2)
% %             plot(bleb_y_filt)
% %             pause
            [peakv,peakh]=min(bleb_y_filt(idxes));
            peakv=-(peakv-bleb_y_baseline);
            peakh=peakh+idxes(1)-1;
            [baselinepeakv,~]=max(bleb_y_filt(baselineidxes));
            baselinepeakv=(baselinepeakv-bleb_y_baseline);
            hafampl=bleb_y_baseline-peakv/2;
            firsth=peakh;
            while firsth>1 & bleb_y_filt(firsth-1)<hafampl
                firsth=firsth-1;
            end
            secondh=peakh;
            while secondh<length(bleb_y_filt) & bleb_y_filt(secondh+1)<hafampl
                secondh=secondh+1;
            end
            bleb_halfwidth=(secondh-firsth)*APdata(api).si;
            maxdeviation=max(abs([max(bleb_y_filt)-bleb_y_baseline,min(bleb_y_filt)-bleb_y_baseline]));
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
        APdata(api).bleb_halfwidth=bleb_halfwidth;
        APdata(api).bleb_SNR=peakv/baselinepeakv;
        APdata(api).maxdeviationontrace=maxdeviation;
    end
    
    baselineamplitudes=sort([APdata(strcmp({APdata.bleb_Amplifiermode},recordingmode)).bleb_baselineamplitude]);
    
    [nbase,xbincenters]=hist(baselineamplitudes,100);
    baselineamplitudes=baselineamplitudes(round(baselinepercentilerange(1)*length(baselineamplitudes)):round(baselinepercentilerange(2)*length(baselineamplitudes)));
    [nbase2,~]=hist(baselineamplitudes,xbincenters);
    minamplitude=mean(baselineamplitudes)+std(baselineamplitudes)*baselinestdmin;
    minamplitude_for_noise=mean(baselineamplitudes)+std(baselineamplitudes)*baselinestdmin_for_noise; 
    todel=find(([APdata.bleb_halfwidth]<minhalfwidth | [APdata.maxdeviationontrace]>minamplitude_for_noise | [APdata.bleb_amplitude]<minamplitude) & strcmp({APdata.bleb_Amplifiermode},recordingmode) );
    
    APdata(todel)=[];
    
%     close all
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
            needed=find(needed&[APdata.axonalAP]==1 & [APdata.amplitude]<.06 & [APdata.halfwidth]>300*10^-6); %
        elseif i==1 % persistent spike
            cim='axonal spike';
            needed=find(needed&[APdata.axonalAP]==1 & [APdata.amplitude]>.07 &  [APdata.baselineval]<-.00); %& [APdata.halfwidth]>300*10^-6
        elseif i==2 % somatic spike
%             cim='somatic spike - hyperpolarized v0';
%             needed=find(needed&[APdata.somaticAP]==1 & [APdata.amplitude]>.07&  [APdata.baselineval]<-.07);
            cim='all axonal';
            needed=find(needed&[APdata.axonalAP]==1);
        elseif i==3 % somatic spike
            cim='somatic spike - depolarized v0';
            needed=find(needed&[APdata.somaticAP]==1 & [APdata.amplitude]>.07&  [APdata.baselineval]>-.05);
        end
        figure(1)
        hsom(i+1)=subplot(4,5,i*5+1);
        hold on
        plot([APdata(needed).time],[APdata(needed).ic_y])
        plot(mean([APdata(needed).time],2),mean([APdata(needed).ic_y],2),'k-','LineWidth',3)
        title(cim)
        hax(i+1)=subplot(4,5,i*5+2);
        hold on
        plot([APdata(needed).time],bsxfun(@(a,b)a-b,[APdata(needed).bleb_y_filt],[APdata(needed).bleb_y_baseline]))
        plot(mean([APdata(needed).time],2),mean(bsxfun(@(a,b)a-b,[APdata(needed).bleb_y_filt],[APdata(needed).bleb_y_baseline]),2),'k-','LineWidth',3)
        if isfield(APdata,'bleb_amplitude')
            [~,amplbins]=hist([APdata(strcmp({APdata.bleb_Amplifiermode},recordingmode)).bleb_amplitude],50);
            [~,ampltbins]=hist([APdata(strcmp({APdata.bleb_Amplifiermode},recordingmode)).bleb_amplitude_t],50);
            [~,hwbins]=hist([APdata(strcmp({APdata.bleb_Amplifiermode},recordingmode)).bleb_halfwidth],50);
            subplot(4,5,i*5+3);
            hist([APdata(needed).bleb_amplitude],amplbins)
            xlabel('peak amplitude')
            subplot(4,5,i*5+4);
            hist([APdata(needed).bleb_amplitude_t],ampltbins)
            xlabel('peak latency')
            subplot(4,5,i*5+5);
            hist([APdata(needed).bleb_halfwidth],hwbins)
            xlabel('bleb halfwidth')
        end
%         linkaxes(hsom,'xy')
%         linkaxes(hax,'xy')
        subplot(4,5,i*5+2);
        axis tight
        subplot(4,5,i*5+1);
        axis tight
        if isfield(APdata,'bleb_amplitude')
            figure(2)
            if i<3
                colorka='rx';
            else
                colorka='ko';
            end
            subplot(5,1,1)
            hold on
            plot([APdata(needed).maxtime],[APdata(needed).bleb_amplitude],colorka)
            ylabel('peak amplitude')
            subplot(5,1,2)
            hold on
            plot([APdata(needed).maxtime],[APdata(needed).bleb_amplitude_t],colorka)
            ylabel('peak latency')
            xlabel('time')
            subplot(5,1,3)
            hold on
            plot([APdata(needed).maxtime],[APdata(needed).bleb_halfwidth],colorka)
            ylabel('halfwidth')
            xlabel('time')
            subplot(5,1,4)
            hold on
            plot([APdata(needed).maxtime],[APdata(needed).bleb_SNR],colorka)
            ylabel('SNR')
            xlabel('time')
            subplot(5,1,5)
            hold on
            plot([APdata(needed).bleb_amplitude_t],[APdata(needed).bleb_amplitude],colorka)
            xlabel('peak latency')
            ylabel('peak amplitude')
        end
    end
    %% run through bleb APs
    window=10;%APnum
    step=5;%APnum
    for i=1:round((length(APdata)-window)/step)
        needed=[1:window]+(i-1)*step;
        figure(3)
        clf
        subplot(3,1,1)
        plot(mean([APdata(needed).time],2),mean([APdata(needed).ic_y],2),'k-','LineWidth',3)
        subplot(3,1,2)
        plot(mean([APdata(needed).time],2),mean(bsxfun(@(a,b)a-b,[APdata(needed).bleb_y_filt],[APdata(needed).bleb_y_baseline]),2),'k-','LineWidth',3)
        subplot(3,1,3)
        hold on
        plot([APdata.maxtime],[APdata.bleb_amplitude],'ko')
        plot([APdata(needed).maxtime],[APdata(needed).bleb_amplitude],'rx','MarkerSize',12)
        pause
    end
    
    
    
end