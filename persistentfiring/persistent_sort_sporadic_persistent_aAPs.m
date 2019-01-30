function eventdata=persistent_sort_sporadic_persistent_aAPs(eventdata,valtozok)
window_def=3;
max_sporadicAP_freq_def=1;
persistent_firing_minaAP_freq_def=10;
valtozok.sAPnumszorzo=2; 
if nargin==1
    window=window_def;
    max_sporadicAP_freq=max_sporadicAP_freq_def;%Hz
    persistent_firing_minaAP_freq=persistent_firing_minaAP_freq_def;%Hz in a window
else
    if isfield(valtozok,'window')
        window=valtozok.window;
    else
        window=window_def;
    end
    if isfield(valtozok,'max_sporadicAP_freq')
        max_sporadicAP_freq=valtozok.max_sporadicAP_freq;
    else
        max_sporadicAP_freq=max_sporadicAP_freq_def;
    end
    if isfield(valtozok,'persistent_firing_minaAP_freq')
        persistent_firing_minaAP_freq=valtozok.persistent_firing_minaAP_freq;
    else
        persistent_firing_minaAP_freq=persistent_firing_minaAP_freq_def;
    end
end
starttime=eventdata(1).sweeptime;
endtime=eventdata(end).maxtime;%
reclength=endtime-starttime;
aAPidxes=find([eventdata.axonalAP]);%&~ [eventdata.stimulated]
aAPtimes=[eventdata(aAPidxes).maxtime];
sAPidxes=find([eventdata.somaticAP] );%&~ [eventdata.stimulated]
sAPtimes=[eventdata(sAPidxes).maxtime];

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
if any([eventdata.aAPfreq_in_window]>persistent_firing_minaAP_freq)
    pfstarttime=eventdata(find([eventdata.aAPfreq_in_window]>persistent_firing_minaAP_freq,1,'first')).maxtime;
else
    pfstarttime=endtime;
end

sporadicaAPidx=[eventdata.axonalAP] & [eventdata.aAPfreq_in_window]>=[eventdata.sAPfreq_in_window]*valtozok.sAPnumszorzo & [eventdata.aAPfreq_in_window]<max_sporadicAP_freq & [eventdata.maxtime]<pfstarttime;
persistentaAPidx=[eventdata.axonalAP] & ([eventdata.aAPfreq_in_window]<[eventdata.sAPfreq_in_window] | [eventdata.aAPfreq_in_window]>max_sporadicAP_freq | [eventdata.maxtime]>=pfstarttime);

for eventi=1:length(eventdata)
    if sporadicaAPidx(eventi)==1
        eventdata(eventi).axonalAP_sporadic=1;
    else
        eventdata(eventi).axonalAP_sporadic=0;
    end
    if persistentaAPidx(eventi)==1
        eventdata(eventi).axonalAP_persistent=1;
    else
        eventdata(eventi).axonalAP_persistent=0;
    end
end