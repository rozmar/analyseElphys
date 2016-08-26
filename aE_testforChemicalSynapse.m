function [tracedata,time]=aE_testforChemicalSynapse(valtozok,dirs,pretraces,posttraces,preevents,postevents)
tracedata=struct;
time=[];
apdiffs=diff([preevents.eventdata.maxtime]);
preevents.eventdata(find(apdiffs==0)+1)=[];
apdiffs=diff([preevents.eventdata.maxtime]);
prediffs=[inf,apdiffs];
postdiffs=[apdiffs,inf];
postpostdiffs=[apdiffs(2:end),inf,inf]; %time after the second spike in the paired pulse protocol should be more than this value
if valtozok.pairedpulseneeded==1
    neededapidx=find(prediffs>valtozok.noAPbeforetheevent & postdiffs<valtozok.pairedpulsedelay+valtozok.pairedpulsejitter/2 & postdiffs>valtozok.pairedpulsedelay-valtozok.pairedpulsejitter/2 & postpostdiffs>valtozok.noAPaftertheevent); %criteria to select presynaptic APs
else
    neededapidx=find(prediffs>valtozok.noAPbeforetheevent & postdiffs>valtozok.noAPaftertheevent); %criteria to select presynaptic APs
end
neededapmaxtimes=[preevents.eventdata(neededapidx).maxtime];
neededapmaxhs=[preevents.eventdata(neededapidx).maxh];
neededapsweepnums=[preevents.eventdata(neededapidx).sweepnum];


for preapnum=1:length(neededapmaxtimes)
    postsweepidx=find(pretraces.bridgeddata(neededapsweepnums(preapnum)).realtime==[posttraces.bridgeddata.realtime] | pretraces.bridgeddata(neededapsweepnums(preapnum)).timertime==[posttraces.bridgeddata.timertime]);
    if ~isempty(postsweepidx) & ~any([postevents.eventdata.maxtime]>neededapmaxtimes(preapnum)-valtozok.baselinelength & [postevents.eventdata.maxtime]<neededapmaxtimes(preapnum)+valtozok.psplength)%
        if isempty(fieldnames(tracedata))
            NEXT=1;
        else
            NEXT=length(tracedata)+1;
        end
        fieldek_bridgeddata=fieldnames(pretraces.bridgeddata);
        fieldek_stimdata={'Amplifiermode','preamplnum'};
        for finum=1:length(fieldek_bridgeddata)
            tracedata(NEXT).(['pre_',fieldek_bridgeddata{finum}])=pretraces.bridgeddata(neededapsweepnums(preapnum)).(fieldek_bridgeddata{finum});
            tracedata(NEXT).(['post_',fieldek_bridgeddata{finum}])=posttraces.bridgeddata(postsweepidx).(fieldek_bridgeddata{finum});
            tracedata(NEXT).apmaxh=neededapmaxhs(preapnum);
        end
        for finum=1:length(fieldek_stimdata)
            tracedata(NEXT).(['pre_',fieldek_stimdata{finum}])=pretraces.stimdata(neededapsweepnums(preapnum)).(fieldek_stimdata{finum});
            tracedata(NEXT).(['post_',fieldek_stimdata{finum}])=posttraces.stimdata(postsweepidx).(fieldek_stimdata{finum});
        end
        tracedata(NEXT).(['pre_Amplifiermode'])=pretraces.stimdata(neededapsweepnums(preapnum)).Amplifiermode;
    end
end
if ~isempty(fieldnames(tracedata))
    [si,ia,ic]=unique(round([tracedata.post_si]*10^6));
    si=si/10^6;
    if length(si)>1
        disp(['si gebasz']) % if not all sweeps are recorded with the same sampling freqency, there might be some problem
        si=max(si);
        for sweepnum=1:length(tracedata) % resampling to lowest sampling rate
            if round((tracedata(sweepnum).pre_si*10^6))/10^6<si
                tracedata(sweepnum).apmaxh=ceil(tracedata(sweepnum).apmaxh*round((tracedata(sweepnum).pre_si*10^6))/(si*10^6));
                tracedata(sweepnum).pre_y=resample(tracedata(sweepnum).pre_y,round((tracedata(sweepnum).pre_si*10^6)),si*10^6);
                tracedata(sweepnum).post_y=resample(tracedata(sweepnum).post_y,round((tracedata(sweepnum).pre_si*10^6)),si*10^6);
                tracedata(sweepnum).post_si=si;
                tracedata(sweepnum).pre_si=si;
            end
            
            
        end
        
        disp(['resampled to ',num2str(si)])
    end
    
    % cutting out needed traces and filtering
    ninq=.5/si;
    [b,a] = butter(valtozok.filterorder,valtozok.cutofffreq/ninq,'low');
    stepback=round(valtozok.baselinelength/si);
    stepforward=round(valtozok.psplength/si);
    time=-stepback*si:si:stepforward*si;
    sweeptodel=[];
    for sweepnum=1:length(tracedata)
        if length(tracedata(sweepnum).pre_y)-stepforward<tracedata(sweepnum).apmaxh
            vegh=length(tracedata(sweepnum).pre_y);
            nanvegen=stepforward-(length(tracedata(sweepnum).pre_y)-tracedata(sweepnum).apmaxh);
        else
            vegh=tracedata(sweepnum).apmaxh+stepforward;
            nanvegen=0;
        end
        if stepback>=tracedata(sweepnum).apmaxh
            kezdeth=1;
            nanelejen=stepback-tracedata(sweepnum).apmaxh+1;
        else
            kezdeth=tracedata(sweepnum).apmaxh-stepback;
            nanelejen=0;
        end
        tempy=[ones(1,500)*tracedata(sweepnum).pre_y(1),tracedata(sweepnum).pre_y,ones(1,500)*tracedata(sweepnum).pre_y(end)];
        tempy=filter(b,a,tempy)';
        tracedata(sweepnum).pre_y=tempy(501:length(tracedata(sweepnum).pre_y)+500);
        tempy=[ones(1,500)*tracedata(sweepnum).post_y(1),tracedata(sweepnum).post_y,ones(1,500)*tracedata(sweepnum).post_y(end)];
        tempy=filter(b,a,tempy)';
        tracedata(sweepnum).post_y=tempy(501:length(tracedata(sweepnum).post_y)+500);
        
        
        tracedata(sweepnum).pre_y=[nan(nanelejen,1);tracedata(sweepnum).pre_y(kezdeth:vegh);nan(nanvegen,1)];
        tracedata(sweepnum).post_y=[nan(nanelejen,1);tracedata(sweepnum).post_y(kezdeth:vegh);nan(nanvegen,1)];
        tracedata(sweepnum).apmaxh=tracedata(sweepnum).apmaxh+nanelejen;
        tracedata(sweepnum).post_y0=nanmean(tracedata(sweepnum).post_y(1:stepback));
        tracedata(sweepnum).pre_y0=nanmean(tracedata(sweepnum).pre_y(1:stepback));
        tracedata(sweepnum).post_y0sd=nanstd(tracedata(sweepnum).post_y(1:stepback));
    end
    
    
    
end