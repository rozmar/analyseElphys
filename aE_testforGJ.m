function [tracedataGJ]=aE_testforGJ(valtozok,dirs,pretraces,posttraces,preevents,postevents)
tracedataGJ=struct;
prehypsweeps=struct;

for sweepnum=1:length(pretraces.stimdata)
    currents=pretraces.stimdata(sweepnum).y(pretraces.stimdata(sweepnum).segmenths);
    currdiffs=[0,diff(currents)];
    currdiffs(end)=[];
    currents(end)=[];
    tdiffs=diff(pretraces.stimdata(sweepnum).segmenths)*pretraces.bridgeddata(sweepnum).si;
    if length(currdiffs)>2
        potentialIDXes=find(currdiffs(1:end-1)<valtozok.gj_mincurrampl&tdiffs(1:end-1)>valtozok.gj_minlinelength&currents(1:end-1)<currents(1));
    else
        potentialIDXes=[];
    end
    for potidxi=1:length(potentialIDXes)
        if currdiffs(potentialIDXes(potidxi))==-currdiffs(potentialIDXes(potidxi)+1);
            
            if isempty(fieldnames(prehypsweeps));
                NEXT=1;
            else
                NEXT=length(prehypsweeps)+1;
            end
            prehypsweeps(NEXT).sweepnum=sweepnum;
            prehypsweeps(NEXT).si=pretraces.bridgeddata(sweepnum).si;
            prehypsweeps(NEXT).current=currdiffs(potentialIDXes(potidxi));
            prehypsweeps(NEXT).length=tdiffs(potentialIDXes(potidxi));
            prehypsweeps(NEXT).starth=pretraces.stimdata(sweepnum).segmenths(potentialIDXes(potidxi));
            prehypsweeps(NEXT).endh=pretraces.stimdata(sweepnum).segmenths(potentialIDXes(potidxi)+1);
            prehypsweeps(NEXT).realtime=pretraces.bridgeddata(sweepnum).realtime;
            prehypsweeps(NEXT).timertime=pretraces.bridgeddata(sweepnum).timertime;
            
        end
    end
end





if ~isempty(fieldnames(prehypsweeps))
    for prehypsweepnum=1:length(prehypsweeps)
        postsweepnums=find(prehypsweeps(prehypsweepnum).realtime==[posttraces.bridgeddata.realtime] | prehypsweeps(prehypsweepnum).timertime==[posttraces.bridgeddata.timertime]) ;
        for postsweepnumi=1:length(postsweepnums) % if more channels are recorded at the same time, multiple channels will be present
            postsweepnum=postsweepnums(postsweepnumi);
            if ~any([postevents.eventdata.maxtime]>prehypsweeps(prehypsweepnum).realtime+prehypsweeps(prehypsweepnum).starth*prehypsweeps(prehypsweepnum).si-valtozok.gj_baselinelength & [postevents.eventdata.maxtime]<prehypsweeps(prehypsweepnum).realtime+prehypsweeps(prehypsweepnum).endh*prehypsweeps(prehypsweepnum).si+valtozok.gj_baselinelength)%
                % csak akkor érdekel minket, hogyha nincs beinjektált
                % áram a másik sejtben
                presweepnum=prehypsweeps(prehypsweepnum).sweepnum;
                if strcmp(pretraces.stimdata(presweepnum).Amplifiermode,'C-Clamp') & strcmp(posttraces.stimdata(postsweepnum).Amplifiermode,'C-Clamp') & any(strfind(pretraces.bridgeddata(presweepnum).channellabel,'Vmon')) & any(strfind(posttraces.bridgeddata(postsweepnum).channellabel,'Vmon')) %both cells must be in c-clamp mode, only Vmons anre considered
                    poststim=posttraces.stimdata(postsweepnum).y;
                    if poststim(prehypsweeps(prehypsweepnum).starth-3) == poststim(prehypsweeps(prehypsweepnum).starth+3)
                        if isempty(fieldnames(tracedataGJ))
                            NEXT=1;
                        else
                            NEXT=length(tracedataGJ)+1;
                        end
                        fieldek=fieldnames(pretraces.bridgeddata);
                        for finum=1:length(fieldek)
                            tracedataGJ(NEXT).(['pre_',fieldek{finum}])=pretraces.bridgeddata(prehypsweeps(prehypsweepnum).sweepnum).(fieldek{finum});
                            tracedataGJ(NEXT).(['post_',fieldek{finum}])=posttraces.bridgeddata(postsweepnum).(fieldek{finum});
                            tracedataGJ(NEXT).endh=prehypsweeps(prehypsweepnum).endh;
                            tracedataGJ(NEXT).starth=prehypsweeps(prehypsweepnum).starth;
                            tracedataGJ(NEXT).length=prehypsweeps(prehypsweepnum).length;
                            tracedataGJ(NEXT).current=round(prehypsweeps(prehypsweepnum).current*10^12);
                        end
                    end
                end
            end
        end
    end
end
if ~isempty(fieldnames(tracedataGJ))
    
    
    
    [si,ia,ic]=unique(round([tracedataGJ.post_si]*10^6));
    si=si/10^6;
    if length(si)>1
        disp(['si gebasz']) % if not all sweeps are recorded with the same sampling freqency, there might be some problem
        si=max(si);
        for sweepnum=1:length(tracedataGJ) % resampling sweeps to lowest sampling rate
            if round((tracedataGJ(sweepnum).pre_si*10^6))/10^6<si
                tracedataGJ(sweepnum).starth=ceil(tracedataGJ(sweepnum).starth/si*tracedataGJ(sweepnum).pre_si);
                tracedataGJ(sweepnum).endh=ceil(tracedataGJ(sweepnum).endh/si*tracedataGJ(sweepnum).pre_si);
                tracedataGJ(sweepnum).pre_y=resample(tracedataGJ(sweepnum).pre_y,round((tracedataGJ(sweepnum).pre_si*10^6)),si*10^6);
                tracedataGJ(sweepnum).post_y=resample(tracedataGJ(sweepnum).post_y,round((tracedataGJ(sweepnum).pre_si*10^6)),si*10^6);
                tracedataGJ(sweepnum).post_si=si;
                tracedataGJ(sweepnum).pre_si=si;
            end
        end
        disp(['resampled to ',num2str(si)])
    end
    
    % cutting out needed traces and filtering
    ninq=.5/si;
    [b,a] = butter(valtozok.filterorder,valtozok.cutofffreq/ninq,'low');
    
    for sweepnum=1:length(tracedataGJ)
        stepback=round(valtozok.gj_baselinelength/si);
        stepforward=round(valtozok.gj_baselinelengthend/si);
        time=-stepback*si:si:(stepforward+tracedataGJ(sweepnum).endh-tracedataGJ(sweepnum).starth)*si;
        
        tracedataGJ(sweepnum).time=time;
        if stepback>=tracedataGJ(sweepnum).starth
            filtered=filter(b,a,[tracedataGJ(sweepnum).pre_y(end:-1:1),tracedataGJ(sweepnum).pre_y]);
            tracedataGJ(sweepnum).pre_y=filtered(end-length(tracedataGJ(sweepnum).pre_y):end)';
            filtered=filter(b,a,[tracedataGJ(sweepnum).post_y(end:-1:1),tracedataGJ(sweepnum).post_y]);
            tracedataGJ(sweepnum).post_y=filtered(end-length(tracedataGJ(sweepnum).post_y):end)';
        else
            tracedataGJ(sweepnum).pre_y=filter(b,a,tracedataGJ(sweepnum).pre_y)';
            tracedataGJ(sweepnum).post_y=filter(b,a,tracedataGJ(sweepnum).post_y)';
        end
        
        
        if length(tracedataGJ(sweepnum).pre_y)-stepforward<tracedataGJ(sweepnum).endh
            vegh=length(tracedataGJ(sweepnum).pre_y);
            nanvegen=stepforward-(length(tracedataGJ(sweepnum).pre_y)-tracedataGJ(sweepnum).endh);
        else
            vegh=tracedataGJ(sweepnum).endh+stepforward;
            nanvegen=0;
        end
        if stepback>=tracedataGJ(sweepnum).starth
            kezdeth=1;
            nanelejen=stepback-tracedataGJ(sweepnum).starth;
        else
            kezdeth=tracedataGJ(sweepnum).starth-stepback;
            nanelejen=0;
        end
        
        tracedataGJ(sweepnum).pre_y=[nan(nanelejen,1);tracedataGJ(sweepnum).pre_y(kezdeth:vegh);nan(nanvegen,1)];
        tracedataGJ(sweepnum).post_y=[nan(nanelejen,1);tracedataGJ(sweepnum).post_y(kezdeth:vegh);nan(nanvegen,1)];
        tracedataGJ(sweepnum).post_y0=nanmean(tracedataGJ(sweepnum).post_y(1:stepback));
        tracedataGJ(sweepnum).pre_y0=nanmean(tracedataGJ(sweepnum).pre_y(1:stepback));
        tracedataGJ(sweepnum).post_y0sd=nanstd(tracedataGJ(sweepnum).post_y(1:stepback));
        %                     figure(333)
        %                     clf
        %                     subplot(2,1,1)
        %                     plot(tracedataGJ(sweepnum).time,tracedataGJ(sweepnum).pre_y);
        %                     subplot(2,1,2)
        %                     plot(tracedataGJ(sweepnum).time,tracedataGJ(sweepnum).post_y);
        %                     pause
    end
    neededtraces=find([tracedataGJ.length]==mode([tracedataGJ.length])&[tracedataGJ.current]==mode([tracedataGJ.current]));
    tracedataGJ=tracedataGJ(neededtraces);
end