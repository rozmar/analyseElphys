function xlsdata = persistent_aAPstatistics(dirs,xlsdata)
brainstatenames={'Slow wave sleep', 'REM sleep', 'Quiet wakefulness','Active wakefulness'};
brainstatenames_varname={'Slow_wave_sleep', 'REM_sleep', 'Quiet_wakefulness','Active_wakefulness'};
% stepsize=1;
window=10;
max_sporadicAP_freq=1;%Hz
persistent_firing_minaAP_freq=2;%Hz in a window
%% aAP and sAP statistics in each experiment
sortedeventdir=[dirs.eventdir,'sorted/'];

files=dir(sortedeventdir);
files([files.isdir])=[];
APdata=struct;
fieldek=fieldnames(xlsdata);

for xlsi=1:length(xlsdata)
    if  xlsdata(xlsi).field==0 & xlsdata(xlsi).juxta==0 %any(strfind(xlsdata(xlsi).anaesthesia,'awake')) &
        if isempty(fieldnames(APdata))
            NEXT=1;
        else
            NEXT=length(APdata)+1;
        end
        for xlsfieldi=1:length(fieldek)
            APdata(NEXT).(fieldek{xlsfieldi})=xlsdata(xlsi).(fieldek{xlsfieldi});
        end
        a=dir([sortedeventdir,xlsdata(xlsi).ID,'.mat']);
        if isempty(a)
            APdata(NEXT).aAPnum=NaN;
            APdata(NEXT).sAPnum=NaN;
            xlsdata(xlsi).axonalAPnum=NaN;
            xlsdata(xlsi).somaticAPnum=NaN;
            xlsdata(xlsi).axonalAPnum_sporadic=NaN;
            xlsdata(xlsi).axonalAPnum_persistent=NaN;
            xlsdata(xlsi).aAPratio=NaN;
            xlsdata(xlsi).recordinglength=NaN;
            xlsdata(xlsi).axonalAPfreq=NaN;
            xlsdata(xlsi).somaticAPfreq=NaN;
            xlsdata(xlsi).aAPisi_mean=NaN;
            xlsdata(xlsi).aAPisi_min=NaN;
            xlsdata(xlsi).aAPisi_max=NaN;
            xlsdata(xlsi).aAPisi_med=NaN;
            xlsdata(xlsi).sAPisi_mean=NaN;
            xlsdata(xlsi).sAPisi_min=NaN;
            xlsdata(xlsi).sAPisi_max=NaN;
            xlsdata(xlsi).sAPisi_med=NaN;
            xlsdata(xlsi).stimulatedsAPnum=NaN;
            for brainstatei=1:length(brainstatenames)
                brainstatename_var=brainstatenames_varname{brainstatei};
                xlsdata(xlsi).([(brainstatename_var),'_timespent'])=NaN;
                xlsdata(xlsi).([(brainstatename_var),'_aAPnum'])=NaN;
                xlsdata(xlsi).([(brainstatename_var),'_aAPfreq'])=NaN;
            end
            xlsdata(xlsi).brainstateless_timespent=NaN;
            xlsdata(xlsi).brainstateless_aAPnum=NaN;
            xlsdata(xlsi).brainstateless_aAPfreq=NaN;
        else
            disp(xlsdata(xlsi).ID)
            load([sortedeventdir,xlsdata(xlsi).ID],'eventdata');
            load([dirs.bridgeddir,xlsdata(xlsi).ID],'bridgeddata');
            %%
            
            starttime=bridgeddata(1).realtime;
            endtime=bridgeddata(end).realtime+bridgeddata(end).si*length(bridgeddata(end).y);
            reclength=endtime-starttime;
            aAPidxes=find([eventdata.axonalAP]);%&~ [eventdata.stimulated]
            aAPtimes=[eventdata(aAPidxes).maxtime];
            sAPidxes=find([eventdata.somaticAP] );%&~ [eventdata.stimulated]
            sAPtimes=[eventdata(sAPidxes).maxtime];
            %%
            for eventi=1:length(eventdata)
                if strcmp(eventdata(eventi).type,'AP')
                    timestart=max(starttime,eventdata(eventi).maxtime-window/2);
                    timeend=min(endtime,eventdata(eventi).maxtime+window/2);
                    timeinterval=timeend-timestart;
                    aAPnum=sum(aAPtimes>timestart & aAPtimes<=timeend);
                    eventdata(eventi).aAPfreq_in_window=aAPnum/timeinterval;
                    sAPnum=sum(sAPtimes>timestart & sAPtimes<=timeend);
                    eventdata(eventi).sAPfreq_in_window=sAPnum/timeinterval;
                else
                    eventdata(eventi).aAPfreq_in_window=NaN;
                    eventdata(eventi).sAPfreq_in_window=NaN;
                end
            end

            %%
            if any([eventdata.aAPfreq_in_window]>persistent_firing_minaAP_freq)
                pfstarttime=eventdata(find([eventdata.aAPfreq_in_window]>persistent_firing_minaAP_freq,1,'first')).maxtime;
            else
                pfstarttime=endtime;
            end
            sporadicaAPidx=[eventdata.axonalAP] & [eventdata.aAPfreq_in_window]>=[eventdata.sAPfreq_in_window] & [eventdata.aAPfreq_in_window]<max_sporadicAP_freq & [eventdata.maxtime]<pfstarttime;
            persistentaAPidx=[eventdata.axonalAP] & ([eventdata.aAPfreq_in_window]<[eventdata.sAPfreq_in_window] | [eventdata.aAPfreq_in_window]>max_sporadicAP_freq | [eventdata.maxtime]>=pfstarttime);
            APdata(NEXT).aAPnum=sum([eventdata.axonalAP]&~[eventdata.stimulated]);
            APdata(NEXT).sAPnum=sum([eventdata.somaticAP]&~[eventdata.stimulated]);
            xlsdata(xlsi).axonalAPnum=APdata(NEXT).aAPnum;
            xlsdata(xlsi).somaticAPnum=APdata(NEXT).sAPnum;
            xlsdata(xlsi).axonalAPnum_sporadic=sum(sporadicaAPidx);
            xlsdata(xlsi).axonalAPnum_persistent=sum(persistentaAPidx);
            xlsdata(xlsi).aAPratio=xlsdata(xlsi).axonalAPnum/(xlsdata(xlsi).axonalAPnum+xlsdata(xlsi).somaticAPnum);
            xlsdata(xlsi).recordinglength=endtime-starttime;
            xlsdata(xlsi).axonalAPfreq=APdata(NEXT).aAPnum/xlsdata(xlsi).recordinglength;
            aAPisis=diff([eventdata([eventdata.axonalAP] & ~[eventdata.stimulated]).maxtime]);
            xlsdata(xlsi).stimulatedsAPnum=sum([[eventdata.somaticAP] & [eventdata.stimulated]]);
            xlsdata(xlsi).somaticAPfreq=APdata(NEXT).sAPnum/xlsdata(xlsi).recordinglength;
            xlsdata(xlsi).aAPtimes=[eventdata([eventdata.axonalAP]&~[eventdata.stimulated]).maxtime]-starttime;
            xlsdata(xlsi).sAPtimes=[eventdata([eventdata.somaticAP]&~[eventdata.stimulated]).maxtime]-starttime;
            if ~isempty(aAPisis)
                xlsdata(xlsi).aAPisi_mean=mean(aAPisis);
                xlsdata(xlsi).aAPisi_min=min(aAPisis);
                xlsdata(xlsi).aAPisi_max=max(aAPisis);
                xlsdata(xlsi).aAPisi_med=median(aAPisis);
            else
                xlsdata(xlsi).aAPisi_mean=NaN;
                xlsdata(xlsi).aAPisi_min=NaN;
                xlsdata(xlsi).aAPisi_max=NaN;
                xlsdata(xlsi).aAPisi_med=NaN;
            end
            sAPisis=diff([eventdata([eventdata.somaticAP] & ~[eventdata.stimulated]).maxtime]);
            if ~isempty(aAPisis)
                xlsdata(xlsi).sAPisi_mean=mean(sAPisis);
                xlsdata(xlsi).sAPisi_min=min(sAPisis);
                xlsdata(xlsi).sAPisi_max=max(sAPisis);
                xlsdata(xlsi).sAPisi_med=median(sAPisis);
            else
                xlsdata(xlsi).sAPisi_mean=NaN;
                xlsdata(xlsi).sAPisi_min=NaN;
                xlsdata(xlsi).sAPisi_max=NaN;
                xlsdata(xlsi).sAPisi_med=NaN;
            end
        end
        if isfield(dirs,'brainstatedir')
        a=dir([dirs.brainstatedir,xlsdata(xlsi).ID,'.mat']);
        if ~isempty(a)
            %%
            sporadicaAPidx=[eventdata.axonalAP] & [eventdata.aAPfreq_in_window]>=[eventdata.sAPfreq_in_window] & [eventdata.aAPfreq_in_window]<max_sporadicAP_freq & [eventdata.maxtime]<pfstarttime;
            persistentaAPidx=[eventdata.axonalAP] & ([eventdata.aAPfreq_in_window]<[eventdata.sAPfreq_in_window] | [eventdata.aAPfreq_in_window]>max_sporadicAP_freq | [eventdata.maxtime]>=pfstarttime);
            load([dirs.brainstatedir,xlsdata(xlsi).ID]);
            aAPpresentsofar=0;
            timespentsofar=0;
            for brainstatei=1:length(brainstatenames)
                brainstatename=brainstatenames{brainstatei};
                brainstatename_var=brainstatenames_varname{brainstatei};
                needed = find(strcmp({BrainStateData.name},brainstatename));
                if ~isempty(needed)
                    timespent=sum([BrainStateData(needed).endtime]-[BrainStateData(needed).starttime]);
%                     aAPtimes=[eventdata([eventdata.axonalAP]&~[eventdata.stimulated]).maxtime];
                    aAPtimes=[eventdata(sporadicaAPidx).maxtime];
                    isAPpresent=false(size(aAPtimes));
                    for brainstateidx=1:length(needed)
                        isAPpresent=isAPpresent| (aAPtimes>[BrainStateData(needed(brainstateidx)).starttime] & aAPtimes<[BrainStateData(needed(brainstateidx)).endtime]);
                    end
                    aAPpresentsofar=aAPpresentsofar+sum(isAPpresent);
                    timespentsofar=timespentsofar+timespent;
                    xlsdata(xlsi).([(brainstatename_var),'_timespent'])=timespent;
                    xlsdata(xlsi).([(brainstatename_var),'_aAPnum'])=sum(isAPpresent);
                    xlsdata(xlsi).([(brainstatename_var),'_aAPfreq'])=sum(isAPpresent)/timespent;
                else
                    xlsdata(xlsi).([(brainstatename_var),'_timespent'])=0;
                    xlsdata(xlsi).([(brainstatename_var),'_aAPnum'])=NaN;
                    xlsdata(xlsi).([(brainstatename_var),'_aAPfreq'])=NaN;
                end
                
            end
            xlsdata(xlsi).brainstateless_timespent=xlsdata(xlsi).recordinglength-timespentsofar;
            xlsdata(xlsi).brainstateless_aAPnum=xlsdata(xlsi).axonalAPnum_sporadic-aAPpresentsofar;
            xlsdata(xlsi).brainstateless_aAPfreq=xlsdata(xlsi).brainstateless_aAPnum/xlsdata(xlsi).brainstateless_timespent;
            %%
%             disp('lol')
        else
            for brainstatei=1:length(brainstatenames)
                brainstatename_var=brainstatenames_varname{brainstatei};
                xlsdata(xlsi).([(brainstatename_var),'_timespent'])=NaN;
                xlsdata(xlsi).([(brainstatename_var),'_aAPnum'])=NaN;
                xlsdata(xlsi).([(brainstatename_var),'_aAPfreq'])=NaN;
            end
            xlsdata(xlsi).brainstateless_timespent=NaN;
            xlsdata(xlsi).brainstateless_aAPnum=NaN;
            xlsdata(xlsi).brainstateless_aAPfreq=NaN;
        end
    end
    else
        xlsdata(xlsi).axonalAPnum=NaN;
        xlsdata(xlsi).somaticAPnum=NaN;
        xlsdata(xlsi).axonalAPnum_sporadic=NaN;
        xlsdata(xlsi).axonalAPnum_persistent=NaN;
        xlsdata(xlsi).aAPratio=NaN;
        xlsdata(xlsi).recordinglength=NaN;
        xlsdata(xlsi).axonalAPfreq=NaN;
        xlsdata(xlsi).somaticAPfreq=NaN;
        xlsdata(xlsi).aAPisi_mean=NaN;
        xlsdata(xlsi).aAPisi_min=NaN;
        xlsdata(xlsi).aAPisi_max=NaN;
        xlsdata(xlsi).aAPisi_med=NaN;
        
        xlsdata(xlsi).sAPisi_mean=NaN;
        xlsdata(xlsi).sAPisi_min=NaN;
        xlsdata(xlsi).sAPisi_max=NaN;
        xlsdata(xlsi).sAPisi_med=NaN;
        xlsdata(xlsi).stimulatedsAPnum=NaN;
        
         for brainstatei=1:length(brainstatenames)
                brainstatename_var=brainstatenames_varname{brainstatei};
                xlsdata(xlsi).([(brainstatename_var),'_timespent'])=NaN;
                xlsdata(xlsi).([(brainstatename_var),'_aAPnum'])=NaN;
                xlsdata(xlsi).([(brainstatename_var),'_aAPfreq'])=NaN;
            end
            xlsdata(xlsi).brainstateless_timespent=NaN;
            xlsdata(xlsi).brainstateless_aAPnum=NaN;
            xlsdata(xlsi).brainstateless_aAPfreq=NaN;
    end
end