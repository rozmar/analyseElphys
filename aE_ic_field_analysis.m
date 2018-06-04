function dataout=aE_ic_field_analysis(dirs,xlsdata,icxlsnum,timeborders,type,additionaldata)
dataout=struct;
plotthestuff=0;
%% field and IC analysis values

% %% slice
% % oscillating slice - sporadic axonal APs
% fieldxlsnum=find(strcmp({xlsdata.ID},'1608041tm_3_1_2'));%341;
% icxlsnum=find(strcmp({xlsdata.ID},'1608041tm_3_5,17_4'));%342;
% timeborders=[0,300]+49821;
% 
% % oscillating slice - rhythmic persistent firing
% fieldxlsnum=find(strcmp({xlsdata.ID},'1611091tm_1_1_3'));%401;
% icxlsnum=find(strcmp({xlsdata.ID},'1611091tm_3_1_4'));%403;
% timeborders=[44800,44900];
% 
% %rat slow oscillation L1 interneuron receiving EPSPs
% fieldxlsnum=find(strcmp({xlsdata.ID},'1611151tm_3_2_3'));
% icxlsnum=find(strcmp({xlsdata.ID},'1611151tm_3_1_4'));
% timeborders=[];
% % 
% %rhythmic persisntent - dsfc - no field
% fieldxlsnum=find(strcmp({xlsdata.ID},'1702101tm_2_33_3'));
% fieldxlsnum=NaN;
% icxlsnum=find(strcmp({xlsdata.ID},'1702101tm_2_1,2_4'));
% timeborders=[53800, 55500];
% 
% %rhythmic persisntent - dsfc - no field - NBQX-GBZ
% fieldxlsnum=NaN;
% icxlsnum=find(strcmp({xlsdata.ID},'1703032tm_2_1_4'));
% timeborders=[59220, 59220+350];
%  %% invivo
% 
%  %nonNGF - no field
% fieldxlsnum=NaN;
% icxlsnum=find(strcmp({xlsdata.ID},'1704191rm_2_1,2_1'));
% timeborders=[67732, 68300];
% 
% %axonal spikes - no field
% fieldxlsnum=NaN;
% icxlsnum=find(strcmp({xlsdata.ID},'1704191rm_3_1,2_1'));
% timeborders=[0, 10000000];
%  
% %  %axonal spikes - no field
% % fieldxlsnum=NaN;
% % icxlsnum=find(strcmp({xlsdata.ID},'1704241rm_1_5,7,9_1'));
% % timeborders=[0,67101010600];
% [icxlsnum,ok] = listdlg('ListString',{xlsdata.ID},'ListSize',[300 600]); % 

% timeborders=[71320 71450];
% timeborders=[0 70700];

%% loading
fieldxlsnum=find(strcmp(xlsdata(icxlsnum).HEKAfname,{xlsdata.HEKAfname}) & [xlsdata.field]==1);
if isempty(fieldxlsnum)
    fieldxlsnum=(NaN);
end

ic=load([dirs.bridgeddir,xlsdata(icxlsnum).ID]);
if isnan(fieldxlsnum)
    field=ic;
else
    field=load([dirs.bridgeddir,xlsdata(fieldxlsnum).ID]);
    jo=0; % -if there is no overlap between the field and the ic recording..
    for i=1:length(ic.lightdata)
        if any([ic.lightdata(i).realtime]==[field.lightdata.realtime])
            jo=1;
        end
    end
    if jo==0
        fieldxlsnum=NaN;
        field=ic;
    end
end
if nargin<6 & ~any(strcmp(fieldnames(additionaldata),'eventdata'))
    load([dirs.eventdir,xlsdata(icxlsnum).ID],'eventdata');
else
    eventdata=additionaldata.eventdata;
end
if ~isempty(timeborders)
    needed=find([eventdata.maxtime]>=timeborders(1)& [eventdata.maxtime]<=timeborders(2));
    eventdata=eventdata(needed);
end
% needed=strcmp({eventdata.type},eventtype);% & ~[eventdata.stimulated] & [eventdata.baselineval]<-.06;
% needed=strcmp({eventdata.type},'ep');% & ~[eventdata.stimulated] & [eventdata.baselineval]<-.06;

% eventdata=eventdata(needed);
diffs=[inf,diff([eventdata.maxtime])];
diffs2=[diff([eventdata.maxtime]),inf];
diffmin=min([diffs;diffs2]);
needed=diffmin>0;
eventdata=eventdata(needed);
apdata=eventdata(strcmp({eventdata.type},'AP') & [eventdata.stimulated]==0);%& [eventdata.baselineval]<-.05
aapdata=eventdata(strcmp({eventdata.type},'AP') & [eventdata.axonalAP] & [eventdata.stimulated]==0);%& [eventdata.baselineval]<-.05
sapdata=eventdata(strcmp({eventdata.type},'AP') & [eventdata.somaticAP] & [eventdata.stimulated]==0);%& [eventdata.baselineval]<-.05
epdata=eventdata(strcmp({eventdata.type},'ep'));
ipdata=eventdata(strcmp({eventdata.type},'ip'));
%
% close all
localextremumwin=.1;
cutofffreq=[1 5];
timebefore=.5;
timeafter=.5;
timebefore_corr=1;
timeafter_corr=1;
NEXT=0;
X=[];
Y=[];
YFIELD=[];
YFIELDreal=[];
Yphase=[];
YfieldAmplitude=[];
YfieldTimetoTrough=[];
YfieldTimetoPeak=[];
APtime=[];
V0=[];
halfperiodlength=[];
ISI=[];
ISI_prev=[];
ISI_next=[];
NEXT=0;
FieldData=struct;
EventTimesRelative=struct;

% analysis relative to the Field
progressbar('analyzing field data')
for fieldsweepnum= 1: length(field.bridgeddata)
    progressbar(fieldsweepnum/ length(field.bridgeddata));
    if isempty(timeborders) | (field.bridgeddata(fieldsweepnum).realtime>=timeborders(1) & field.bridgeddata(fieldsweepnum).realtime<=timeborders(2))
    si=field.bridgeddata(fieldsweepnum).si;
    [b,a]=butter(1,cutofffreq/(1/field.bridgeddata(fieldsweepnum).si)/2,'bandpass');
    [bb,aa]=butter(1,500/(1/field.bridgeddata(fieldsweepnum).si)/2,'low');
    [bic,aic]=butter(1,1500/(1/field.bridgeddata(fieldsweepnum).si)/2,'low');
    yfield=filtfilt(b,a,field.bridgeddata(fieldsweepnum).y);
    
    yfieldorig=yfield;
    yfieldtoshow=filtfilt(bb,aa,field.bridgeddata(fieldsweepnum).y);
    stepback=round(timebefore/si);
    stepforward=round(timeafter/si);
    sweepnum=find([ic.bridgeddata.realtime]==field.bridgeddata(fieldsweepnum).realtime);
    if ~isempty(sweepnum)
        
        ic.bridgeddata(sweepnum).yfield=field.bridgeddata(fieldsweepnum).y;
        ic.bridgeddata(sweepnum).yfield_filt=yfield;
    end
    lestep=round(localextremumwin/si);
    minh2=stepback;
    maxh=stepback;
    maxh2=stepback;
    if ~isempty(sweepnum)%% analysis relative to the Field
        if isnan(fieldxlsnum) | strcmp(type,'peak') % the waveform is upside-down if the recording is intracellular
            yfield=yfield*-1;
        end
        icy=ic.bridgeddata(sweepnum).y;
        sweeptime=[0:length(icy)-1]*si+ic.bridgeddata(sweepnum).realtime;
        while maxh2<length(yfield)-stepforward-stepback
            maxh=maxh2;
            while maxh<length(yfield)-lestep & max(yfield(maxh:maxh+lestep))>max(yfield(maxh-lestep:maxh))
                maxh=maxh+lestep;
            end
            if maxh+2*lestep<length(yfield)
                minh2=maxh+lestep;
            end
            while minh2<length(yfield)-2*lestep & minh2<length(yfield)-2*stepforward& min(yfield(minh2:minh2+lestep))<min(yfield(minh2-lestep:minh2))
                minh2=minh2+lestep;
            end
            if minh2+2*lestep<length(yfield)
                maxh2=minh2+lestep;
            end
            while maxh2<length(yfield)-2*lestep & max(yfield(maxh2:maxh2+lestep))>max(yfield(maxh2-lestep:maxh2))
                maxh2=maxh2+lestep;
            end
            
            [maxval,sech]=max(yfield(maxh-lestep:maxh));
            sech=sech+maxh-lestep;
            [minval,firsth2]=min(yfield(minh2-lestep:minh2));
            firsth2=firsth2+minh2-lestep;
            [maxval2,sech2]=max(yfield(maxh2:maxh2+lestep));
            sech2=sech2+maxh2;
            
            if isnan(fieldxlsnum) | strcmp(type,'peak') % the waveform is upside-down if the recording is intracellular
                maxval=maxval*-1;
                minval=minval*-1;
                maxval2=maxval2*-1;
            end
            
            time=[-stepback:stepforward]*si;
            fieldidx=[-stepback:stepforward]+firsth2;
            if isnan(fieldxlsnum) | strcmp(type,'peak') % the waveform is upside-down if the recording is intracellular
                fieldy=-1*yfield(fieldidx);
            else
                fieldy=yfield(fieldidx);
            end
            fieldtoshow=yfieldtoshow(fieldidx);
            intracell=icy(fieldidx);
            intracell=filtfilt(bic,aic,intracell);
            relativepeakhs=[sech,firsth2,sech2]-firsth2+stepback;
            relativepeakhs(relativepeakhs<1)=1;
            relativepeakhs(relativepeakhs>stepback+stepforward)=stepback+stepforward;
            NEXT=NEXT+1;
            amplitudes=abs(diff(fieldy(relativepeakhs)));
            FieldData(NEXT).time=time';
            FieldData(NEXT).si=mode(diff(time));
            FieldData(NEXT).fieldtodetect=fieldy';
            FieldData(NEXT).fieldfilt=fieldtoshow';
            FieldData(NEXT).ic=intracell';
            FieldData(NEXT).amplitudes=amplitudes;
            FieldData(NEXT).meanamplitude=nanmean(amplitudes);
            FieldData(NEXT).maxamplitude=nanmax(amplitudes);
            FieldData(NEXT).peakvalue1=fieldy(relativepeakhs(1));
            FieldData(NEXT).troughvalue=fieldy(relativepeakhs(2));
            FieldData(NEXT).peakvalue2=fieldy(relativepeakhs(3));
            FieldData(NEXT).relativepeakhs=relativepeakhs;
            FieldData(NEXT).prevpeak=sech;
            FieldData(NEXT).troughh=firsth2;
            FieldData(NEXT).nextpeak=sech2;
            FieldData(NEXT).fieldsweepnum=fieldsweepnum;
            FieldData(NEXT).sweepnum=sweepnum;
            FieldData(NEXT).troughtime=sweeptime(firsth2);
            
            eventdata=apdata;
            neededevents=find([eventdata.maxtime]>=FieldData(NEXT).troughtime-timebefore & [eventdata.maxtime]<=FieldData(NEXT).troughtime+timeafter);
            FieldData(NEXT).relativeAPtimes=[eventdata(neededevents).maxtime]-FieldData(NEXT).troughtime;
            if isempty(FieldData(NEXT).relativeAPtimes)
                FieldData(NEXT).apnum=0;
            else
                FieldData(NEXT).apnum=length(FieldData(NEXT).relativeAPtimes);
            end
            
            eventdata=aapdata;
            neededevents=find([eventdata.maxtime]>=FieldData(NEXT).troughtime-timebefore & [eventdata.maxtime]<=FieldData(NEXT).troughtime+timeafter);
            FieldData(NEXT).relativeaxonalAPtimes=[eventdata(neededevents).maxtime]-FieldData(NEXT).troughtime;
            if isempty(FieldData(NEXT).relativeaxonalAPtimes)
                FieldData(NEXT).axonalapnum=0;
            else
                FieldData(NEXT).axonalapnum=length(FieldData(NEXT).relativeaxonalAPtimes);
            end
            
            eventdata=epdata;
            neededevents=find([eventdata.maxtime]>=FieldData(NEXT).troughtime-timebefore & [eventdata.maxtime]<=FieldData(NEXT).troughtime+timeafter);
            FieldData(NEXT).relativeeptimes=[eventdata(neededevents).maxtime]-FieldData(NEXT).troughtime;
            if isempty(FieldData(NEXT).relativeeptimes)
                FieldData(NEXT).epnum=0;
            else
                FieldData(NEXT).epnum=length(FieldData(NEXT).relativeeptimes);
            end
            
            eventdata=ipdata;
            neededevents=find([eventdata.maxtime]>=FieldData(NEXT).troughtime-timebefore & [eventdata.maxtime]<=FieldData(NEXT).troughtime+timeafter);
            FieldData(NEXT).relativeiptimes=[eventdata(neededevents).maxtime]-FieldData(NEXT).troughtime;
            if isempty(FieldData(NEXT).relativeiptimes)
                FieldData(NEXT).ipnum=0;
            else
                FieldData(NEXT).ipnum=length(FieldData(NEXT).relativeiptimes);
            end
            FieldData(NEXT).eventnum=FieldData(NEXT).ipnum+FieldData(NEXT).epnum+FieldData(NEXT).apnum;
            
%             figure(3)
%             clf
%             subplot(3,1,1)
%             plot(time,fieldy,'k-')
%             hold on;
%             plot(time(relativepeakhs),fieldy(relativepeakhs),'ro')
%             axis tight
%             subplot(3,1,2)
%             plot(time,fieldtoshow,'k-')
%             hold on;
%             plot(time(relativepeakhs),fieldtoshow(relativepeakhs),'ro')
%             axis tight
%             subplot(3,1,3)
%             plot(time,intracell,'k-')
%             axis tight
%             pause
        end
%         if NEXT>0
%         figure(2)
%         clf
%         needed=[FieldData.fieldsweepnum]==fieldsweepnum;
%          subplot(3,1,1)
%             plot(sweeptime,yfieldorig,'k-')
%             hold on;
%             plot(sweeptime([FieldData(needed).troughh]),yfieldorig([FieldData(needed).troughh]),'ro')
%             axis tight
%             subplot(3,1,2)
%             plot(sweeptime,yfieldtoshow,'k-')
%             hold on;
%             plot(sweeptime([FieldData(needed).troughh]),yfieldtoshow([FieldData(needed).troughh]),'ro')
%             axis tight
%             axis tight
%             subplot(3,1,3)
%             plot(sweeptime,icy,'k-')
%             axis tight
% %             pause
%         end
    end
    end
end
FieldDataoriginal=FieldData;
%% calculating cross correlations
% events.AP=apdata;
% events.aAP=aapdata;
% events.sAP=sapdata;
% events.EP=epdata;
% events.IP=ipdata;
% tocalculate=fieldnames(events);
% for i=1:length(tocalculate)
%     events_now=events.(tocalculate{i});
%     for j=1:length(tocalculate)
%         events_now2=events.(tocalculate{j});
%         for eventi=1:length(events_now);
%             maxtime=events_now(eventi).maxtime;
%             neededevents=find([events_now2.maxtime]>=maxtime-timebefore_corr & [events_now2.maxtime]<=maxtime+timeafter_corr);
%             data=[events_now2(neededevents).maxtime]-maxtime;
%             data(data==0)=[];
%             events_now(eventi).(['relativetimes_',tocalculate{j}])=data;
%         end
%     end
%     events.(tocalculate{i})=events_now;
% end
% 
% dataout.events=events;
dataout.FieldData=FieldData;
dataout.bridgeddata=ic.bridgeddata;
%% plotting
FieldData=FieldDataoriginal;
for i=1:length(FieldData)
    FieldData(i).medicV=median(FieldData(i).ic);
end
FieldData([FieldData.eventnum]==0)=[];
timestep=.05;
bins=[-timebefore:timestep:timeafter];
bins_corr=[-timebefore_corr:timestep:timeafter_corr];
timewindowforfieldamplitude=10;
% needed=[FieldData.maxamplitude]>=2*10^-5&[FieldData.maxamplitude]<=10*10^-5;
% FieldData=FieldData(needed);
for i = 1:length(FieldData)
    needed=FieldData(i).troughtime-timewindowforfieldamplitude/2<=[FieldData.troughtime] & FieldData(i).troughtime+timewindowforfieldamplitude/2>=[FieldData.troughtime];
    FieldData(i).medianamplitude=median([FieldData(needed).maxamplitude]);
end
dataout.FieldData=FieldData;
needed=[FieldData.maxamplitude]>=[FieldData.medianamplitude]/3;%&[FieldData.maxamplitude]<=[FieldData.medianamplitude]*5;% & [FieldData.apnum]>0 ;
FieldData=FieldData(needed);
periodlength=nanmean([NaN,diff([FieldData.troughtime]);diff([FieldData.troughtime,NaN])]);
needed=periodlength>.2 & periodlength<2;
FieldData=FieldData(needed);
% needed=[FieldData.medicV]<-.055;
% FieldData=FieldData(needed);
% ido=[FieldData.time];%bsxfun(@(A,B) A+B,[FieldData.time],[FieldData.troughtime]);
%%
if plotthestuff==1
    figure(1)
    clf
    plot([FieldData.troughtime],[FieldData.maxamplitude],'k-')
    hold on
    % plot([FieldData.troughtime],[FieldData.medianamplitude],'r-')
    figure(2)
    clf
    sis=unique([FieldData.si]);
    for i=1:length(sis)
        needed=[FieldData.si]==sis(i);
        h(1)=subplot(4,1,1);
        plot([FieldData(needed).time],[FieldData(needed).ic]*1000)
        hold on
        % plot(ido,mean([FieldData.ic],2)*1000,'k-','LineWIdth',3);
        ylabel('IC (mV)')
        axis tight
        h(2)=subplot(4,1,2);
        hold on
        plot([FieldData(needed).time],[FieldData(needed).fieldtodetect]*1000); %bsxfun(@(A,B )A-B,[FieldData.fieldtodetect],[FieldData.troughvalue])
        plot([FieldData(needed).time],mean([FieldData(needed).fieldtodetect],2)*1000,'k-','LineWIdth',3);
        if isnan(fieldxlsnum)
            ylabel('Filtered intra (mV) -no field ')
        else
            ylabel('Field (mV)')
        end
        axis tight
    end
    h(3)=subplot(4,1,3);
    [nAP,bins]=hist([FieldData.relativeAPtimes],[-timebefore:timestep:timeafter]);
    [naxonalAP,~]=hist([FieldData.relativeaxonalAPtimes],[-timebefore:timestep:timeafter]);
    hold on
    bar(bins,nAP,'b')
    bar(bins,naxonalAP,'r')
    ylabel('AP count')
    axis tight
    h(4)=subplot(4,1,4);
    hist([FieldData.relativeeptimes],[-timebefore:timestep:timeafter])
    axis tight
    ylabel('EPSP count')
    xlabel('time (s)')
    linkaxes(h,'x');
    %%
    figure(3)
    clf
    for i=1:length(tocalculate)
        for j=1:length(tocalculate)
            subplot(length(tocalculate),length(tocalculate),(i-1)*length(tocalculate)+j);
            if length((events.(tocalculate{i})))>1
                hist([events.(tocalculate{i}).(['relativetimes_',tocalculate{j}])],bins_corr);
                axis tight
            end
            
            title([tocalculate{j},' relative to ',tocalculate{i}])
        end
    end
    
    %% saving images
    
    ID=xlsdata(icxlsnum).ID;
    figure(2)
    saveas(gcf,[dirs.figuresdir,ID,'_field_anal'],'pdf')
    saveas(gcf,[dirs.figuresdir,ID,'_field_anal'],'jpg')
    figure(3)
    saveas(gcf,[dirs.figuresdir,ID,'_event_correlations'],'pdf')
    saveas(gcf,[dirs.figuresdir,ID,'_event_correlations'],'jpg')
end
disp('field analysis finished')
%%
return
% %%
% events.AP=apdata;
% events.ep=epdata;
% events.ip=ipdata;
% eventtypes=fieldnames(events);
% for eventtypei=1:1%length(eventtypes)
%     eventdata=events.(eventtypes{eventtypei});
%     prevsweepnum=0;
%     for eventi=1:length(eventdata)
%         si=eventdata(eventi).si;
%         stepback=round(timebefore/si);
%         stepforward=round(timeafter/si);
%         maxh=eventdata(eventi).maxh;
%         sweepnum=eventdata(eventi).sweepnum;
%         lestep=round(localextremumwin/si);
%         if maxh>stepback+lestep & length(ic.bridgeddata(sweepnum).y)>maxh+stepforward
%             x=[-stepback:stepforward]*si;
%             if sweepnum~=prevsweepnum
%                 yicraw=ic.bridgeddata(sweepnum).y;
%                 fieldsweepnum=find([field.bridgeddata.realtime]==eventdata(eventi).sweeptime);
%                 yfieldraw=field.bridgeddata(fieldsweepnum).y;
%                 [b_deltaFilt,a_deltaFilt] = ellip(1,0.1,40,[0.1 4].*si);
%                 yfield_deltafilt = filtfilt(b_deltaFilt,a_deltaFilt,yfieldraw);
%                 yfield_hil = hilbert(yfield_deltafilt);
%                 K_yfield = abs(yfield_hil);
%                 yfield_deltafilt_phi= angle(yfield_hil); 
%                 prevsweepnum=sweepnum;
%             end  
%             yicraw_now=yicraw(maxh-stepback:maxh+stepforward);
%             yfieldraw_now=yfieldraw(maxh-stepback:maxh+stepforward);
%             yfield_deltafilt_now=yfield_deltafilt(maxh-stepback:maxh+stepforward);
%             maxtime=eventdata(eventi).maxtime;
%             
%             %%
%             pause
% %             NINCSENKÉSZ!!!!! asdfasfgasfsa
%             prevtroughidx=find([FieldData.troughtime]<maxtime,1,'last');
%             prevtroughtime=FieldData(prevtroughidx).troughtime;
%             nexttroughidx=find([FieldData.troughtime]>maxtime,1,'first');
%             nexttroughtime=FieldData(nexttroughidx).troughtime;
%             
%             sweeptimevector=[0:length(ic.bridgeddata(sweepnum).y)]*ic.bridgeddata(sweepnum).si+ic.bridgeddata(sweepnum).realtime;
%             peaktimes=sort(sweeptimevector([FieldData(prevtroughidx).prevpeak,FieldData(prevtroughidx).nextpeak,FieldData(nexttroughidx).prevpeak,FieldData(nexttroughidx).nextpeak]));
%             
%             prevpeaktime=peaktimes(find(peaktimes<maxtime,1,'last'));
%             nextpeaktime=peaktimes(find(peaktimes>maxtime,1,'first'));
%             
%             timetoprevpeak=prevpeaktime-maxtime;
%             timetonextpeak=nextpeaktime-maxtime;
%             timetoprevtrough=prevtroughtime-maxtime;
%             timetonexttrough=nexttroughtime-maxtime;
%             
%             
%             X(NEXT,:)=x;
%             if x(1)==0
%                 disp('lol')
%             end
%             Y(NEXT,:)=y;
%             YFIELD(NEXT,:)=yfield;
%             YFIELDreal(NEXT,:)=yfieldorig*1000000;
%             %         Yphase(NEXT)=phase;
%             YfieldAmplitude(NEXT)=fieldamplitude;
%             YfieldTimetoPrevTrough(NEXT)=timetoprevtrough;
%             YfieldTimetoNextTrough(NEXT)=timetonexttrough;
%             %         YfieldTimetoPeak(NEXT)=timetopeak;
%             APtime(NEXT)=eventdata(eventi).maxtime;
%             V0(NEXT)=eventdata(eventi).baselineval;
%             ISI(NEXT)=diffmin(eventi);
%             ISI_prev(NEXT)=diffs(eventi);
%             ISI_next(NEXT)=diffs2(eventi);
%         else
%             disp('lol')
%         end
%         progressbar(eventi/length(eventdata));
%     end
% end
% 
% %%
% return
% %% analysis relative to the APs
% eventdata=apdata;
% for eventi=1:length(eventdata)
%     si=eventdata(eventi).si;
%     stepback=round(timebefore/si);
%     stepforward=round(timeafter/si);
%     maxh=eventdata(eventi).maxh;
%     sweepnum=eventdata(eventi).sweepnum;
%     lestep=round(localextremumwin/si);
%     if maxh>stepback+lestep & length(ic.bridgeddata(sweepnum).y)>maxh+stepforward
%         NEXT=NEXT+1;
%         x=[-stepback:stepforward]*si;
%         %         y=nan(size(x));
%         y=ic.bridgeddata(sweepnum).y(maxh-stepback:maxh+stepforward);
%         fieldsweepnum=find([field.bridgeddata.realtime]==eventdata(eventi).sweeptime);
%         [b,a]=butter(1,cutofffreq/(1/field.bridgeddata(fieldsweepnum).si)/2,'bandpass');
%         [bb,aa]=butter(1,500/(1/field.bridgeddata(fieldsweepnum).si)/2,'low');
%         
%         yfield=filtfilt(b,a,field.bridgeddata(fieldsweepnum).y);
%         yfieldorig=yfield;
%         yfieldorig=yfieldorig(maxh-stepback:maxh+stepforward);
%         yfield=yfield(maxh-stepback:maxh+stepforward);
%         yfieldtoshow=filtfilt(bb,aa,field.bridgeddata(fieldsweepnum).y);
%         yfieldtoshow=yfieldtoshow(maxh-stepback:maxh+stepforward);
%         
%         
%         maxtime=eventdata(eventi).maxtime;
%         prevtroughidx=find([FieldData.troughtime]<maxtime,1,'last');
%         prevtroughtime=FieldData(prevtroughidx).troughtime;
%         nexttroughidx=find([FieldData.troughtime]>maxtime,1,'first');
%         nexttroughtime=FieldData(nexttroughidx).troughtime;
%         
%         sweeptimevector=[0:length(ic.bridgeddata(sweepnum).y)]*ic.bridgeddata(sweepnum).si+ic.bridgeddata(sweepnum).realtime;
%         peaktimes=sort(sweeptimevector([FieldData(prevtroughidx).prevpeak,FieldData(prevtroughidx).nextpeak,FieldData(nexttroughidx).prevpeak,FieldData(nexttroughidx).nextpeak]));
%         
%         prevpeaktime=peaktimes(find(peaktimes<maxtime,1,'last'));
%         nextpeaktime=peaktimes(find(peaktimes>maxtime,1,'first'));
%         
%         timetoprevpeak=prevpeaktime-maxtime;
%         timetonextpeak=nextpeaktime-maxtime;
%         timetoprevtrough=prevtroughtime-maxtime;
%         timetonexttrough=nexttroughtime-maxtime;
%         
% 
%         X(NEXT,:)=x;
%         if x(1)==0
%             disp('lol')
%         end
%         Y(NEXT,:)=y;
%         YFIELD(NEXT,:)=yfield;
%         YFIELDreal(NEXT,:)=yfieldorig*1000000;
% %         Yphase(NEXT)=phase;
%         YfieldAmplitude(NEXT)=fieldamplitude;
%         YfieldTimetoPrevTrough(NEXT)=timetoprevtrough;
%         YfieldTimetoNextTrough(NEXT)=timetonexttrough;
% %         YfieldTimetoPeak(NEXT)=timetopeak;
%         APtime(NEXT)=eventdata(eventi).maxtime;
%         V0(NEXT)=eventdata(eventi).baselineval;
%         ISI(NEXT)=diffmin(eventi);
%         ISI_prev(NEXT)=diffs(eventi);
%         ISI_next(NEXT)=diffs2(eventi);
%     else
%         disp('lol')
%     end
%     progressbar(eventi/length(eventdata));
% end
% %%
% periodlength=YfieldTimetoNextTrough-YfieldTimetoPrevTrough;
% % periodlength=halfperiodlength;
% minamplitude=0E-5;
% minperiodlength=.1;%         peaktimes=sort([prevtroughtime+FieldData(prevtroughidx).time(FieldData(prevtroughidx).prevpeak),prevtroughtime+FieldData(prevtroughidx).time(FieldData(prevtroughidx).nextpeak),nexttroughtime+FieldData(nexttroughidx).time(FieldData(nexttroughidx).prevpeak),nexttroughtime+FieldData(nexttroughidx).time(FieldData(nexttroughidx).nextpeak)]);
% 
% limit=.3;
% lfplimit=.5;
% close all
% figure(1)
% clf
% subplot(5,4,1)
% needed=ISI>limit & YfieldAmplitude>minamplitude & periodlength>minperiodlength;
% rose(degtorad(Yphase(needed)),20);
% [tout,rout] =rose(degtorad(Yphase(needed)),20);
% figure(2)
% plotCircularBarchart(rout, 0, 'solitary APs')
% figure(1)
% title('solitary APs')
% subplot(5,4,5)
% % shadedErrorBar(mean(X(needed,:)),mean(YFIELDreal(needed,:)),std(YFIELDreal(needed,:)),'k')
% hold on
% plot(mean(X(needed,:)),mean(YFIELDreal(needed,:)),'k-','LineWidth',3)
% axis tight
% % ylim([-1 1])
% % xlim([-lfplimit lfplimit])
% subplot(5,4,9)
% hold all
% plot(X(needed,:)',Y(needed,:)')
% xlim([-limit limit])
% subplot(5,4,13);
% prevtroughtimes=[YfieldTimetoPrevTrough(needed)];
% nexttroughtimes=[YfieldTimetoNextTrough(needed)];
% periodlengths=nexttroughtimes-prevtroughtimes;
% periodprigresses=-prevtroughtimes./periodlengths;
% hist(-[prevtroughtimes./periodlengths,nexttroughtimes./periodlengths],[-1:.1:1])
% xlabel('normalized time to trough')
% subplot(5,4,17);
% plot(APtime(needed),periodprigresses,'ko-')
% xlabel('Time (s)')
% ylabel('Phase (trought-to-trough')
% 
% subplot(5,4,2)
% needed=ISI<limit  & YfieldAmplitude>minamplitude& periodlength>minperiodlength ;
% rose(degtorad(Yphase(needed)),20);
% [tout,rout] =rose(degtorad(Yphase(needed)),20);
% figure(3)
% plotCircularBarchart(rout, 0, 'grouped APs')
% figure(1)
% title('grouped APs')
% subplot(5,4,6)
% % shadedErrorBar(mean(X(needed,:)),mean(YFIELDreal(needed,:)),std(YFIELDreal(needed,:)),'k')
% hold on
% plot(mean(X(needed,:)),mean(YFIELDreal(needed,:)),'k-','LineWidth',3)
% axis tight
% % ylim([-1 1])
% % xlim([-lfplimit lfplimit])
% subplot(5,4,10)
% hold all
% plot(X(needed,:)',Y(needed,:)')
% xlim([-limit limit])
% subplot(5,4,14);
% prevtroughtimes=[YfieldTimetoPrevTrough(needed)];
% nexttroughtimes=[YfieldTimetoNextTrough(needed)];
% periodlengths=nexttroughtimes-prevtroughtimes;
% periodprigresses=-prevtroughtimes./periodlengths;
% hist(-[prevtroughtimes./periodlengths,nexttroughtimes./periodlengths],[-1:.1:1])
% xlabel('normalized time to trough')
% subplot(5,4,18);
% plot(APtime(needed),periodprigresses,'ko-')
% xlabel('Time (s)')
% ylabel('Phase (trought-to-trough')
% 
% subplot(5,4,3)
% needed=ISI_prev>limit & ISI_next<limit & YfieldAmplitude>minamplitude& periodlength>minperiodlength ;
% rose(degtorad(Yphase(needed)),20);
% [tout,rout] =rose(degtorad(Yphase(needed)),20);
% figure(4)
% plotCircularBarchart(rout, 0, 'group leading APs')
% figure(1)
% title('group leading APs')% %%
% events.AP=apdata;
% events.ep=epdata;
% events.ip=ipdata;
% eventtypes=fieldnames(events);
% for eventtypei=1:1%length(eventtypes)
%     eventdata=events.(eventtypes{eventtypei});
%     prevsweepnum=0;
%     for eventi=1:length(eventdata)
%         si=eventdata(eventi).si;
%         stepback=round(timebefore/si);
%         stepforward=round(timeafter/si);
%         maxh=eventdata(eventi).maxh;
%         sweepnum=eventdata(eventi).sweepnum;
%         lestep=round(localextremumwin/si);
%         if maxh>stepback+lestep & length(ic.bridgeddata(sweepnum).y)>maxh+stepforward
%             x=[-stepback:stepforward]*si;
%             if sweepnum~=prevsweepnum
%                 yicraw=ic.bridgeddata(sweepnum).y;
%                 fieldsweepnum=find([field.bridgeddata.realtime]==eventdata(eventi).sweeptime);
%                 yfieldraw=field.bridgeddata(fieldsweepnum).y;
%                 [b_deltaFilt,a_deltaFilt] = ellip(1,0.1,40,[0.1 4].*si);
%                 yfield_deltafilt = filtfilt(b_deltaFilt,a_deltaFilt,yfieldraw);
%                 yfield_hil = hilbert(yfield_deltafilt);
%                 K_yfield = abs(yfield_hil);
%                 yfield_deltafilt_phi= angle(yfield_hil); 
%                 prevsweepnum=sweepnum;
%             end  
%             yicraw_now=yicraw(maxh-stepback:maxh+stepforward);
%             yfieldraw_now=yfieldraw(maxh-stepback:maxh+stepforward);
%             yfield_deltafilt_now=yfield_deltafilt(maxh-stepback:maxh+stepforward);
%             maxtime=eventdata(eventi).maxtime;
%             
%             %%
%             pause
% %             NINCSENKÉSZ!!!!! asdfasfgasfsa
%             prevtroughidx=find([FieldData.troughtime]<maxtime,1,'last');
%             prevtroughtime=FieldData(prevtroughidx).troughtime;
%             nexttroughidx=find([FieldData.troughtime]>maxtime,1,'first');
%             nexttroughtime=FieldData(nexttroughidx).troughtime;
%             
%             sweeptimevector=[0:length(ic.bridgeddata(sweepnum).y)]*ic.bridgeddata(sweepnum).si+ic.bridgeddata(sweepnum).realtime;
%             peaktimes=sort(sweeptimevector([FieldData(prevtroughidx).prevpeak,FieldData(prevtroughidx).nextpeak,FieldData(nexttroughidx).prevpeak,FieldData(nexttroughidx).nextpeak]));
%             
%             prevpeaktime=peaktimes(find(peaktimes<maxtime,1,'last'));
%             nextpeaktime=peaktimes(find(peaktimes>maxtime,1,'first'));
%             
%             timetoprevpeak=prevpeaktime-maxtime;
%             timetonextpeak=nextpeaktime-maxtime;
%             timetoprevtrough=prevtroughtime-maxtime;
%             timetonexttrough=nexttroughtime-maxtime;
%             
%             
%             X(NEXT,:)=x;
%             if x(1)==0
%                 disp('lol')
%             end
%             Y(NEXT,:)=y;
%             YFIELD(NEXT,:)=yfield;
%             YFIELDreal(NEXT,:)=yfieldorig*1000000;
%             %         Yphase(NEXT)=phase;
%             YfieldAmplitude(NEXT)=fieldamplitude;
%             YfieldTimetoPrevTrough(NEXT)=timetoprevtrough;
%             YfieldTimetoNextTrough(NEXT)=timetonexttrough;
%             %         YfieldTimetoPeak(NEXT)=timetopeak;
%             APtime(NEXT)=eventdata(eventi).maxtime;
%             V0(NEXT)=eventdata(eventi).baselineval;
%             ISI(NEXT)=diffmin(eventi);
%             ISI_prev(NEXT)=diffs(eventi);
%             ISI_next(NEXT)=diffs2(eventi);
%         else
%             disp('lol')
%         end
%         progressbar(eventi/length(eventdata));
%     end
% end
% 
% %%
% return
% %% analysis relative to the APs
% eventdata=apdata;
% for eventi=1:length(eventdata)
%     si=eventdata(eventi).si;
%     stepback=round(timebefore/si);
%     stepforward=round(timeafter/si);
%     maxh=eventdata(eventi).maxh;
%     sweepnum=eventdata(eventi).sweepnum;
%     lestep=round(localextremumwin/si);
%     if maxh>stepback+lestep & length(ic.bridgeddata(sweepnum).y)>maxh+stepforward
%         NEXT=NEXT+1;
%         x=[-stepback:stepforward]*si;
%         %         y=nan(size(x));
%         y=ic.bridgeddata(sweepnum).y(maxh-stepback:maxh+stepforward);
%         fieldsweepnum=find([field.bridgeddata.realtime]==eventdata(eventi).sweeptime);
%         [b,a]=butter(1,cutofffreq/(1/field.bridgeddata(fieldsweepnum).si)/2,'bandpass');
%         [bb,aa]=butter(1,500/(1/field.bridgeddata(fieldsweepnum).si)/2,'low');
%         
%         yfield=filtfilt(b,a,field.bridgeddata(fieldsweepnum).y);
%         yfieldorig=yfield;
%         yfieldorig=yfieldorig(maxh-stepback:maxh+stepforward);
%         yfield=yfield(maxh-stepback:maxh+stepforward);
%         yfieldtoshow=filtfilt(bb,aa,field.bridgeddata(fieldsweepnum).y);
%         yfieldtoshow=yfieldtoshow(maxh-stepback:maxh+stepforward);
%         
%         
%         maxtime=eventdata(eventi).maxtime;
%         prevtroughidx=find([FieldData.troughtime]<maxtime,1,'last');
%         prevtroughtime=FieldData(prevtroughidx).troughtime;
%         nexttroughidx=find([FieldData.troughtime]>maxtime,1,'first');
%         nexttroughtime=FieldData(nexttroughidx).troughtime;
%         
%         sweeptimevector=[0:length(ic.bridgeddata(sweepnum).y)]*ic.bridgeddata(sweepnum).si+ic.bridgeddata(sweepnum).realtime;
%         peaktimes=sort(sweeptimevector([FieldData(prevtroughidx).prevpeak,FieldData(prevtroughidx).nextpeak,FieldData(nexttroughidx).prevpeak,FieldData(nexttroughidx).nextpeak]));
%         
%         prevpeaktime=peaktimes(find(peaktimes<maxtime,1,'last'));
%         nextpeaktime=peaktimes(find(peaktimes>maxtime,1,'first'));
%         
%         timetoprevpeak=prevpeaktime-maxtime;
%         timetonextpeak=nextpeaktime-maxtime;
%         timetoprevtrough=prevtroughtime-maxtime;
%         timetonexttrough=nexttroughtime-maxtime;
%         
% 
%         X(NEXT,:)=x;
%         if x(1)==0
%             disp('lol')
%         end
%         Y(NEXT,:)=y;
%         YFIELD(NEXT,:)=yfield;
%         YFIELDreal(NEXT,:)=yfieldorig*1000000;
% %         Yphase(NEXT)=phase;
%         YfieldAmplitude(NEXT)=fieldamplitude;
%         YfieldTimetoPrevTrough(NEXT)=timetoprevtrough;
%         YfieldTimetoNextTrough(NEXT)=timetonexttrough;
% %         YfieldTimetoPeak(NEXT)=timetopeak;
%         APtime(NEXT)=eventdata(eventi).maxtime;
%         V0(NEXT)=eventdata(eventi).baselineval;
%         ISI(NEXT)=diffmin(eventi);
%         ISI_prev(NEXT)=diffs(eventi);
%         ISI_next(NEXT)=diffs2(eventi);
%     else
%         disp('lol')
%     end
%     progressbar(eventi/length(eventdata));
% end
% %%
% periodlength=YfieldTimetoNextTrough-YfieldTimetoPrevTrough;
% % periodlength=halfperiodlength;
% minamplitude=0E-5;
% minperiodlength=.1;%         peaktimes=sort([prevtroughtime+FieldData(prevtroughidx).time(FieldData(prevtroughidx).prevpeak),prevtroughtime+FieldData(prevtroughidx).time(FieldData(prevtroughidx).nextpeak),nexttroughtime+FieldData(nexttroughidx).time(FieldData(nexttroughidx).prevpeak),nexttroughtime+FieldData(nexttroughidx).time(FieldData(nexttroughidx).nextpeak)]);
% 
% limit=.3;
% lfplimit=.5;
% close all
% figure(1)
% clf
% subplot(5,4,1)
% needed=ISI>limit & YfieldAmplitude>minamplitude & periodlength>minperiodlength;
% rose(degtorad(Yphase(needed)),20);
% [tout,rout] =rose(degtorad(Yphase(needed)),20);
% figure(2)
% plotCircularBarchart(rout, 0, 'solitary APs')
% figure(1)
% title('solitary APs')
% subplot(5,4,5)
% % shadedErrorBar(mean(X(needed,:)),mean(YFIELDreal(needed,:)),std(YFIELDreal(needed,:)),'k')
% hold on
% plot(mean(X(needed,:)),mean(YFIELDreal(needed,:)),'k-','LineWidth',3)
% axis tight
% % ylim([-1 1])
% % xlim([-lfplimit lfplimit])
% subplot(5,4,9)
% hold all
% plot(X(needed,:)',Y(needed,:)')
% xlim([-limit limit])
% subplot(5,4,13);
% prevtroughtimes=[YfieldTimetoPrevTrough(needed)];
% nexttroughtimes=[YfieldTimetoNextTrough(needed)];
% periodlengths=nexttroughtimes-prevtroughtimes;
% periodprigresses=-prevtroughtimes./periodlengths;
% hist(-[prevtroughtimes./periodlengths,nexttroughtimes./periodlengths],[-1:.1:1])
% xlabel('normalized time to trough')
% subplot(5,4,17);
% plot(APtime(needed),periodprigresses,'ko-')
% xlabel('Time (s)')
% ylabel('Phase (trought-to-trough')
% 
% subplot(5,4,2)
% needed=ISI<limit  & YfieldAmplitude>minamplitude& periodlength>minperiodlength ;
% rose(degtorad(Yphase(needed)),20);
% [tout,rout] =rose(degtorad(Yphase(needed)),20);
% figure(3)
% plotCircularBarchart(rout, 0, 'grouped APs')
% figure(1)
% title('grouped APs')
% subplot(5,4,6)
% % shadedErrorBar(mean(X(needed,:)),mean(YFIELDreal(needed,:)),std(YFIELDreal(needed,:)),'k')
% hold on
% plot(mean(X(needed,:)),mean(YFIELDreal(needed,:)),'k-','LineWidth',3)
% axis tight
% % ylim([-1 1])
% % xlim([-lfplimit lfplimit])
% subplot(5,4,10)
% hold all
% plot(X(needed,:)',Y(needed,:)')
% xlim([-limit limit])
% subplot(5,4,14);
% prevtroughtimes=[YfieldTimetoPrevTrough(needed)];
% nexttroughtimes=[YfieldTimetoNextTrough(needed)];
% periodlengths=nexttroughtimes-prevtroughtimes;
% periodprigresses=-prevtroughtimes./periodlengths;
% hist(-[prevtroughtimes./periodlengths,nexttroughtimes./periodlengths],[-1:.1:1])
% xlabel('normalized time to trough')
% subplot(5,4,18);
% plot(APtime(needed),periodprigresses,'ko-')
% xlabel('Time (s)')
% ylabel('Phase (trought-to-trough')
% 
% subplot(5,4,3)
% needed=ISI_prev>limit & ISI_next<limit & YfieldAmplitude>minamplitude& periodlength>minperiodlength ;
% rose(degtorad(Yphase(needed)),20);
% [tout,rout] =rose(degtorad(Yphase(needed)),20);
% figure(4)
% plotCircularBarchart(rout, 0, 'group leading APs')
% figure(1)
% title('group leading APs')
% subplot(5,4,7)
% % shadedErrorBar(mean(X(needed,:)),mean(YFIELDreal(needed,:)),std(YFIELDreal(needed,:)),'k')
% hold on
% plot(mean(X(needed,:)),mean(YFIELDreal(needed,:)),'k-','LineWidth',3)
% axis tight
% % ylim([-1 1])persistent_experimental.m
% % xlim([-lfplimit lfplimit])
% subplot(5,4,11)
% hold all
% plot(X(needed,:)',Y(needed,:)')
% xlim([-limit limit])
% subplot(5,4,15);
% prevtroughtimes=[YfieldTimetoPrevTrough(needed)];
% nexttroughtimes=[YfieldTimetoNextTrough(needed)];
% periodlengths=nexttroughtimes-prevtroughtimes;
% periodprigresses=-prevtroughtimes./periodlengths;
% hist(-[prevtroughtimes./periodlengths,nexttroughtimes./periodlengths],[-1:.1:1])
% xlabel('normalized time to trough')
% subplot(5,4,19);
% plot(APtime(needed),periodprigresses,'ko-')
% xlabel('Time (s)')
% ylabel('Phase (trought-to-trough')
% 
% subplot(5,4,4)
% needed=YfieldAmplitude>minamplitude& periodlength>minperiodlength ;
% rose(degtorad(Yphase(needed)),20);
% [tout,rout] =rose(degtorad(Yphase(needed)),20);
% figure(5)
% plotCircularBarchart(rout, 0, 'ALL APs')
% figure(1)
% title('ALL APs')
% subplot(5,4,8)
% % shadedErrorBar(mean(X(needed,:)),mean(YFIELDreal(needed,:)),std(YFIELDreal(needed,:)),'k')
% hold on
% plot(mean(X(needed,:)),mean(YFIELDreal(needed,:)),'k-','LineWidth',3)
% axis tight
% % ylim([-1 1])
% % xlim([-lfplimit lfplimit])
% subplot(5,4,12)
% hold all
% plot(X(needed,:)',Y(needed,:)')
% xlim([-limit limit])
% subplot(5,4,16);
% prevtroughtimes=[YfieldTimetoPrevTrough(needed)];
% nexttroughtimes=[YfieldTimetoNextTrough(needed)];
% periodlengths=nexttroughtimes-prevtroughtimes;
% periodprigresses=-prevtroughtimes./periodlengths;
% hist(-[prevtroughtimes./periodlengths,nexttroughtimes./periodlengths],[-1:.1:1])
% xlabel('normalized time to trough')
% subplot(5,4,20);
% plot(APtime(needed),periodprigresses,'ko-')
% xlabel('Time (s)')
% ylabel('Phase (trought-to-trough')
% % %%
% % close all
% % needed=ISI_prev>limit & ISI_next<limit & YfieldAmplitude>minamplitude& periodlength>minperiodlength ;
% % neededidxs=find(needed);
% % for i=1:length(neededidxs)
% %     idxnow=neededidxs(i);
% %     figure(23)
% %     clf
% %     subplot(3,1,1)
% %     plot(X(idxnow,:),Y(idxnow,:),'k-')
% %     subplot(3,1,2)
% %     plot(X(idxnow,:),YFIELD(idxnow,:),'k-')
% %     axis tight
% %     ylim([-1.5 1.5])
% %     subplot(3,1,3)
% %     plot(X(idxnow,:),YFIELDreal(idxnow,:),'k-')
% %     axis tight
% %     pause% %%
% events.AP=apdata;
% events.ep=epdata;
% events.ip=ipdata;
% eventtypes=fieldnames(events);
% for eventtypei=1:1%length(eventtypes)
%     eventdata=events.(eventtypes{eventtypei});
%     prevsweepnum=0;
%     for eventi=1:length(eventdata)
%         si=eventdata(eventi).si;
%         stepback=round(timebefore/si);
%         stepforward=round(timeafter/si);
%         maxh=eventdata(eventi).maxh;
%         sweepnum=eventdata(eventi).sweepnum;
%         lestep=round(localextremumwin/si);
%         if maxh>stepback+lestep & length(ic.bridgeddata(sweepnum).y)>maxh+stepforward
%             x=[-stepback:stepforward]*si;
%             if sweepnum~=prevsweepnum
%                 yicraw=ic.bridgeddata(sweepnum).y;
%                 fieldsweepnum=find([field.bridgeddata.realtime]==eventdata(eventi).sweeptime);
%                 yfieldraw=field.bridgeddata(fieldsweepnum).y;
%                 [b_deltaFilt,a_deltaFilt] = ellip(1,0.1,40,[0.1 4].*si);
%                 yfield_deltafilt = filtfilt(b_deltaFilt,a_deltaFilt,yfieldraw);
%                 yfield_hil = hilbert(yfield_deltafilt);
%                 K_yfield = abs(yfield_hil);
%                 yfield_deltafilt_phi= angle(yfield_hil); 
%                 prevsweepnum=sweepnum;
%             end  
%             yicraw_now=yicraw(maxh-stepback:maxh+stepforward);
%             yfieldraw_now=yfieldraw(maxh-stepback:maxh+stepforward);
%             yfield_deltafilt_now=yfield_deltafilt(maxh-stepback:maxh+stepforward);
%             maxtime=eventdata(eventi).maxtime;
%             
%             %%
%             pause
% %             NINCSENKÉSZ!!!!! asdfasfgasfsa
%             prevtroughidx=find([FieldData.troughtime]<maxtime,1,'last');
%             prevtroughtime=FieldData(prevtroughidx).troughtime;
%             nexttroughidx=find([FieldData.troughtime]>maxtime,1,'first');
%             nexttroughtime=FieldData(nexttroughidx).troughtime;
%             
%             sweeptimevector=[0:length(ic.bridgeddata(sweepnum).y)]*ic.bridgeddata(sweepnum).si+ic.bridgeddata(sweepnum).realtime;
%             peaktimes=sort(sweeptimevector([FieldData(prevtroughidx).prevpeak,FieldData(prevtroughidx).nextpeak,FieldData(nexttroughidx).prevpeak,FieldData(nexttroughidx).nextpeak]));
%             
%             prevpeaktime=peaktimes(find(peaktimes<maxtime,1,'last'));
%             nextpeaktime=peaktimes(find(peaktimes>maxtime,1,'first'));
%             
%             timetoprevpeak=prevpeaktime-maxtime;
%             timetonextpeak=nextpeaktime-maxtime;
%             timetoprevtrough=prevtroughtime-maxtime;
%             timetonexttrough=nexttroughtime-maxtime;
%             
%             
%             X(NEXT,:)=x;
%             if x(1)==0
%                 disp('lol')
%             end
%             Y(NEXT,:)=y;
%             YFIELD(NEXT,:)=yfield;
%             YFIELDreal(NEXT,:)=yfieldorig*1000000;
%             %         Yphase(NEXT)=phase;
%             YfieldAmplitude(NEXT)=fieldamplitude;
%             YfieldTimetoPrevTrough(NEXT)=timetoprevtrough;
%             YfieldTimetoNextTrough(NEXT)=timetonexttrough;
%             %         YfieldTimetoPeak(NEXT)=timetopeak;
%             APtime(NEXT)=eventdata(eventi).maxtime;
%             V0(NEXT)=eventdata(eventi).baselineval;
%             ISI(NEXT)=diffmin(eventi);
%             ISI_prev(NEXT)=diffs(eventi);
%             ISI_next(NEXT)=diffs2(eventi);
%         else
%             disp('lol')
%         end
%         progressbar(eventi/length(eventdata));
%     end
% end
% 
% %%
% return
% %% analysis relative to the APs
% eventdata=apdata;
% for eventi=1:length(eventdata)
%     si=eventdata(eventi).si;
%     stepback=round(timebefore/si);
%     stepforward=round(timeafter/si);
%     maxh=eventdata(eventi).maxh;
%     sweepnum=eventdata(eventi).sweepnum;
%     lestep=round(localextremumwin/si);
%     if maxh>stepback+lestep & length(ic.bridgeddata(sweepnum).y)>maxh+stepforward
%         NEXT=NEXT+1;
%         x=[-stepback:stepforward]*si;
%         %         y=nan(size(x));
%         y=ic.bridgeddata(sweepnum).y(maxh-stepback:maxh+stepforward);
%         fieldsweepnum=find([field.bridgeddata.realtime]==eventdata(eventi).sweeptime);
%         [b,a]=butter(1,cutofffreq/(1/field.bridgeddata(fieldsweepnum).si)/2,'bandpass');
%         [bb,aa]=butter(1,500/(1/field.bridgeddata(fieldsweepnum).si)/2,'low');
%         
%         yfield=filtfilt(b,a,field.bridgeddata(fieldsweepnum).y);
%         yfieldorig=yfield;
%         yfieldorig=yfieldorig(maxh-stepback:maxh+stepforward);
%         yfield=yfield(maxh-stepback:maxh+stepforward);
%         yfieldtoshow=filtfilt(bb,aa,field.bridgeddata(fieldsweepnum).y);
%         yfieldtoshow=yfieldtoshow(maxh-stepback:maxh+stepforward);
%         
%         
%         maxtime=eventdata(eventi).maxtime;
%         prevtroughidx=find([FieldData.troughtime]<maxtime,1,'last');
%         prevtroughtime=FieldData(prevtroughidx).troughtime;
%         nexttroughidx=find([FieldData.troughtime]>maxtime,1,'first');
%         nexttroughtime=FieldData(nexttroughidx).troughtime;
%         
%         sweeptimevector=[0:length(ic.bridgeddata(sweepnum).y)]*ic.bridgeddata(sweepnum).si+ic.bridgeddata(sweepnum).realtime;
%         peaktimes=sort(sweeptimevector([FieldData(prevtroughidx).prevpeak,FieldData(prevtroughidx).nextpeak,FieldData(nexttroughidx).prevpeak,FieldData(nexttroughidx).nextpeak]));
%         
%         prevpeaktime=peaktimes(find(peaktimes<maxtime,1,'last'));
%         nextpeaktime=peaktimes(find(peaktimes>maxtime,1,'first'));
%         55130
%         timetoprevpeak=prevpeaktime-maxtime;
%         timetonextpeak=nextpeaktime-maxtime;
%         timetoprevtrough=prevtroughtime-maxtime;
%         timetonexttrough=nexttroughtime-maxtime;
%         
% 
%         X(NEXT,:)=x;
%         if x(1)==0
%             disp('lol')
%         end
%         Y(NEXT,:)=y;
%         YFIELD(NEXT,:)=yfield;
%         YFIELDreal(NEXT,:)=yfieldorig*1000000;
% %         Yphase(NEXT)=phase;
%         YfieldAmplitude(NEXT)=fieldamplitude;
%         YfieldTimetoPrevTrough(NEXT)=timetoprevtrough;
%         YfieldTimetoNextTrough(NEXT)=timetonexttrough;
% %         YfieldTimetoPeak(NEXT)=timetopeak;
%         APtime(NEXT)=eventdata(eventi).maxtime;
%         V0(NEXT)=eventdata(eventi).baselineval;
%         ISI(NEXT)=diffmin(eventi);
%         ISI_prev(NEXT)=diffs(eventi);
%         ISI_next(NEXT)=diffs2(eventi);
%     else
%         disp('lol')
%     end
%     progressbar(eventi/length(eventdata));
% end
% %%
% periodlength=YfieldTimetoNextTrough-YfieldTimetoPrevTrough;
% % periodlength=halfperiodlength;
% minamplitude=0E-5;
% minperiodlength=.1;%         peaktimes=sort([prevtroughtime+FieldData(prevtroughidx).time(FieldData(prevtroughidx).prevpeak),prevtroughtime+FieldData(prevtroughidx).time(FieldData(prevtroughidx).nextpeak),nexttroughtime+FieldData(nexttroughidx).time(FieldData(nexttroughidx).prevpeak),nexttroughtime+FieldData(nexttroughidx).time(FieldData(nexttroughidx).nextpeak)]);
% 
% limit=.3;
% lfplimit=.5;
% close all
% figure(1)
% clf
% subplot(5,4,1)
% needed=ISI>limit & YfieldAmplitude>minamplitude & periodlength>minperiodlength;
% rose(degtorad(Yphase(needed)),20);
% [tout,rout] =rose(degtorad(Yphase(needed)),20);
% figure(2)
% plotCircularBarchart(rout, 0, 'solitary APs')
% figure(1)
% title('solitary APs')
% subplot(5,4,5)
% % shadedErrorBar(mean(X(needed,:)),mean(YFIELDreal(needed,:)),std(YFIELDreal(needed,:)),'k')
% hold on
% plot(mean(X(needed,:)),mean(YFIELDreal(needed,:)),'k-','LineWidth',3)
% axis tight
% % ylim([-1 1])
% % xlim([-lfplimit lfplimit])
% subplot(5,4,9)
% hold all
% plot(X(needed,:)',Y(needed,:)')
% xlim([-limit limit])
% subplot(5,4,13);
% prevtroughtimes=[YfieldTimetoPrevTrough(needed)];
% nexttroughtimes=[YfieldTimetoNextTrough(needed)];
% periodlengths=nexttroughtimes-prevtroughtimes;
% periodprigresses=-prevtroughtimes./periodlengths;
% hist(-[prevtroughtimes./periodlengths,nexttroughtimes./periodlengths],[-1:.1:1])
% xlabel('normalized time to trough')
% subplot(5,4,17);
% plot(APtime(needed),periodprigresses,'ko-')
% xlabel('Time (s)')
% ylabel('Phase (trought-to-trough')
% 
% subplot(5,4,2)
% needed=ISI<limit  & YfieldAmplitude>minamplitude& periodlength>minperiodlength ;
% rose(degtorad(Yphase(needed)),20);
% [tout,rout] =rose(degtorad(Yphase(needed)),20);
% figure(3)
% plotCircularBarchart(rout, 0, 'grouped APs')
% figure(1)
% title('grouped APs')
% subplot(5,4,6)
% % shadedErrorBar(mean(X(needed,:)),mean(YFIELDreal(needed,:)),std(YFIELDreal(needed,:)),'k')
% hold on
% plot(mean(X(needed,:)),mean(YFIELDreal(needed,:)),'k-','LineWidth',3)
% axis tight
% % ylim([-1 1])
% % xlim([-lfplimit lfplimit])
% subplot(5,4,10)
% hold all
% plot(X(needed,:)',Y(needed,:)')
% xlim([-limit limit])
% subplot(5,4,14);
% prevtroughtimes=[YfieldTimetoPrevTrough(needed)];
% nexttroughtimes=[YfieldTimetoNextTrough(needed)];
% periodlengths=nexttroughtimes-prevtroughtimes;
% periodprigresses=-prevtroughtimes./periodlengths;
% hist(-[prevtroughtimes./periodlengths,nexttroughtimes./periodlengths],[-1:.1:1])
% xlabel('normalized time to trough')
% subplot(5,4,18);
% plot(APtime(needed),periodprigresses,'ko-')
% xlabel('Time (s)')
% ylabel('Phase (trought-to-trough')
% 
% subplot(5,4,3)
% needed=ISI_prev>limit & ISI_next<limit & YfieldAmplitude>minamplitude& periodlength>minperiodlength ;
% rose(degtorad(Yphase(needed)),20);
% [tout,rout] =rose(degtorad(Yphase(needed)),20);
% figure(4)
% plotCircularBarchart(rout, 0, 'group leading APs')
% figure(1)
% title('group leading APs')
% subplot(5,4,7)
% % shadedErrorBar(mean(X(needed,:)),mean(YFIELDreal(needed,:)),std(YFIELDreal(needed,:)),'k')
% hold on
% plot(mean(X(needed,:)),mean(YFIELDreal(needed,:)),'k-','LineWidth',3)
% axis tight
% % ylim([-1 1])persistent_experimental.m
% % xlim([-lfplimit lfplimit])
% subplot(5,4,11)
% hold all
% plot(X(needed,:)',Y(needed,:)')
% xlim([-limit limit])
% subplot(5,4,15);
% prevtroughtimes=[YfieldTimetoPrevTrough(needed)];
% nexttroughtimes=[YfieldTimetoNextTrough(needed)];
% periodlengths=nexttroughtimes-prevtroughtimes;
% periodprigresses=-prevtroughtimes./periodlengths;
% hist(-[prevtroughtimes./periodlengths,nexttroughtimes./periodlengths],[-1:.1:1])
% xlabel('normalized time to trough')
% subplot(5,4,19);
% plot(APtime(needed),periodprigresses,'ko-')
% xlabel('Time (s)')
% ylabel('Phase (trought-to-trough')
% 
% subplot(5,4,4)
% needed=YfieldAmplitude>minamplitude& periodlength>minperiodlength ;
% rose(degtorad(Yphase(needed)),20);
% [tout,rout] =rose(degtorad(Yphase(needed)),20);
% figure(5)
% plotCircularBarchart(rout, 0, 'ALL APs')
% figure(1)
% title('ALL APs')
% subplot(5,4,8)
% % shadedErrorBar(mean(X(needed,:)),mean(YFIELDreal(needed,:)),std(YFIELDreal(needed,:)),'k')
% hold on
% plot(mean(X(needed,:)),mean(YFIELDreal(needed,:)),'k-','LineWidth',3)
% axis tight
% % ylim([-1 1])
% % xlim([-lfplimit lfplimit])
% subplot(5,4,12)
% hold all
% plot(X(needed,:)',Y(needed,:)')
% xlim([-limit limit])
% subplot(5,4,16);
% prevtroughtimes=[YfieldTimetoPrevTrough(needed)];
% nexttroughtimes=[YfieldTimetoNextTrough(needed)];
% periodlengths=nexttroughtimes-prevtroughtimes;
% periodprigresses=-prevtroughtimes./periodlengths;
% hist(-[prevtroughtimes./periodlengths,nexttroughtimes./periodlengths],[-1:.1:1])
% xlabel('normalized time to trough')
% subplot(5,4,20);
% plot(APtime(needed),periodprigresses,'ko-')
% xlabel('Time (s)')
% ylabel('Phase (trought-to-trough')
% % %%
% % close all
% % needed=ISI_prev>limit & ISI_next<limit & YfieldAmplitude>minamplitude& periodlength>minperiodlength ;
% % neededidxs=find(needed);
% % for i=1:length(neededidxs)
% %     idxnow=neededidxs(i);
% %     figure(23)
% %     clf
% %     subplot(3,1,1)
% %     plot(X(idxnow,:),Y(idxnow,:),'k-')
% %     subplot(3,1,2)
% %     plot(X(idxnow,:),YFIELD(idxnow,:),'k-')
% %     axis tight
% %     ylim([-1.5 1.5])
% %     subplot(3,1,3)
% %     plot(X(idxnow,:),YFIELDreal(idxnow,:),'k-')
% %     axis tight
% %     pause


% subplot(5,4,7)
% % shadedErrorBar(mean(X(needed,:)),mean(YFIELDreal(needed,:)),std(YFIELDreal(needed,:)),'k')
% hold on
% plot(mean(X(needed,:)),mean(YFIELDreal(needed,:)),'k-','LineWidth',3)
% axis tight
% % ylim([-1 1])persistent_experimental.m
% % xlim([-lfplimit lfplimit])
% subplot(5,4,11)
% hold all
% plot(X(needed,:)',Y(needed,:)')
% xlim([-limit limit])
% subplot(5,4,15);
% prevtroughtimes=[YfieldTimetoPrevTrough(needed)];
% nexttroughtimes=[YfieldTimetoNextTrough(needed)];
% periodlengths=nexttroughtimes-prevtroughtimes;
% periodprigresses=-prevtroughtimes./periodlengths;
% hist(-[prevtroughtimes./periodlengths,nexttroughtimes./periodlengths],[-1:.1:1])
% xlabel('normalized time to trough')
% subplot(5,4,19);
% plot(APtime(needed),periodprigresses,'ko-')
% xlabel('Time (s)')
% ylabel('Phase (trought-to-trough')
% 
% subplot(5,4,4)
% needed=YfieldAmplitude>minamplitude& periodlength>minperiodlength ;
% rose(degtorad(Yphase(needed)),20);
% [tout,rout] =rose(degtorad(Yphase(needed)),20);
% figure(5)
% plotCircularBarchart(rout, 0, 'ALL APs')
% figure(1)
% title('ALL APs')
% subplot(5,4,8)
% % shadedErrorBar(mean(X(needed,:)),mean(YFIELDreal(needed,:)),std(YFIELDreal(needed,:)),'k')
% hold on
% plot(mean(X(needed,:)),mean(YFIELDreal(needed,:)),'k-','LineWidth',3)
% axis tight
% % ylim([-1 1])
% % xlim([-lfplimit lfplimit])
% subplot(5,4,12)
% hold all
% plot(X(needed,:)',Y(needed,:)')
% xlim([-limit limit])
% subplot(5,4,16);
% prevtroughtimes=[YfieldTimetoPrevTrough(needed)];
% nexttroughtimes=[YfieldTimetoNextTrough(needed)];
% periodlengths=nexttroughtimes-prevtroughtimes;
% periodprigresses=-prevtroughtimes./periodlengths;
% hist(-[prevtroughtimes./periodlengths,nexttroughtimes./periodlengths],[-1:.1:1])
% xlabel('normalized time to trough')
% subplot(5,4,20);
% plot(APtime(needed),periodprigresses,'ko-')
% xlabel('Time (s)')
% ylabel('Phase (trought-to-trough')
% % %%
% % close all
% % needed=ISI_prev>limit & ISI_next<limit & YfieldAmplitude>minamplitude& periodlength>minperiodlength ;
% % neededidxs=find(needed);
% % for i=1:length(neededidxs)
% %     idxnow=neededidxs(i);
% %     figure(23)
% %     clf
% %     subplot(3,1,1)
% %     plot(X(idxnow,:),Y(idxnow,:),'k-')
% %     subplot(3,1,2)
% %     plot(X(idxnow,:),YFIELD(idxnow,:),'k-')
% %     axis tight
% %     ylim([-1.5 1.5])
% %     subplot(3,1,3)
% %     plot(X(idxnow,:),YFIELDreal(idxnow,:),'k-')
% %     axis tight
% %     pause
% end