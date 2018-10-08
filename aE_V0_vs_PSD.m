function aE_V0_vs_PSD(dirs,xlsdata,xlsnum,valtozok)

if nargin<4
    valtozok=struct;
    valtozok.movingvindowsize=3;
    valtozok.movingvindowstep=.5; %seconds for downsampling and median filtering
    valtozok.timeborders=[0 inf];
    valtozok.frequencyrange=[1 4];
    valtozok.PSDonfield=true;
    valtozok.minsweeplength=valtozok.movingvindowsize/2;
end
% % xlsnum=find(strcmp({xlsdata.ID},'1705301rm_1_2_1'));





fieldxlsnum=find(strcmp(xlsdata(xlsnum).HEKAfname,{xlsdata.HEKAfname}) & [xlsdata.field]==1);
if isempty(fieldxlsnum)
    fieldxlsnum=(NaN);
end

load([dirs.bridgeddir,xlsdata(xlsnum).ID]);

if isnan(fieldxlsnum) | ~valtozok.PSDonfield
    load([dirs.PSDdir,xlsdata(xlsnum).ID]);
    PSDonfield=false;
else
    load([dirs.PSDdir,xlsdata(fieldxlsnum).ID]);
    PSDonfield=true;
end
if isfield(dirs,'brainstatedir')
    a=dir([dirs.brainstatedir,xlsdata(xlsnum).ID,'.mat']);
    if ~isempty(a)
        load([dirs.brainstatedir,xlsdata(xlsnum).ID,'.mat'])
        statenames=unique({BrainStateData.name});
    else
        BrainStateData=[];
        statenames=[];
    end
else
    BrainStateData=[];
        statenames=[];
end

v0_PSD_data=struct;

for sweepi=1:length(bridgeddata)
    hossz=length(bridgeddata(sweepi).y);
    si=bridgeddata(sweepi).si;
    sweepstart=bridgeddata(sweepi).realtime;
    sweepend=bridgeddata(sweepi).realtime+hossz*si;
%     fieldsweepnum=find([PSDdata.realtime]==bridgeddata(sweepi).realtime);
    if si*hossz>valtozok.minsweeplength & (isempty(valtozok.timeborders) | (sweepstart>=valtozok.timeborders(1) & sweepstart<=valtozok.timeborders(2)) | (sweepend>=valtozok.timeborders(1) & sweepend<=valtozok.timeborders(2)) | (sweepstart<=valtozok.timeborders(1) & sweepend>=valtozok.timeborders(2)))% &~isempty(fieldsweepnum)
        %%
        movingvindowstep=round(valtozok.movingvindowstep/si);
        movingvindowsize=round(valtozok.movingvindowsize/si);
        icV=bridgeddata(sweepi).y;
        time=[1:length(icV)]*si-si+bridgeddata(sweepi).realtime;
        icV=medfilt1([icV(ceil(movingvindowsize/2):-1:1),icV,icV(end:-1:end-ceil(movingvindowsize/2))],movingvindowsize);
        icV=icV([1:hossz]+ceil(movingvindowsize/2)+1);
        icV=downsample(icV,movingvindowstep);
        time=downsample(time,movingvindowstep);
        
        fieldsweepnum=find([PSDdata.realtime]==bridgeddata(sweepi).realtime);
        si_PSD=PSDdata(fieldsweepnum).si_powerMatrix;
        frequencyVector=PSDdata(fieldsweepnum).frequencyVector;
        powerMatrix=PSDdata(fieldsweepnum).powerMatrix;
        neededfreqidx=frequencyVector>=min(valtozok.frequencyrange) & frequencyVector<=max(valtozok.frequencyrange);
        powerMatrix=powerMatrix(neededfreqidx,:);
        PSDtime=[1:length(PSDdata(fieldsweepnum).y)]*si_PSD+PSDdata(fieldsweepnum).realtime-si_PSD;
        PSD_max=zeros(size(icV));
        PSD_mean=zeros(size(icV));
        PSD_median=zeros(size(icV));
        stateidx=zeros(size(icV));
        for timei=1:length(time)
            needed=PSDtime>=time(timei)-valtozok.movingvindowsize/2 & PSDtime<=time(timei)+valtozok.movingvindowsize/2;
            datanow=powerMatrix(:,needed);
            PSD_max(timei)=max(datanow(:));
            PSD_mean(timei)=mean(datanow(:));
            PSD_median(timei)=median(datanow(:));
            
            if ~isempty(BrainStateData)
                idx=find(time(timei)>[BrainStateData.starttime] & time(timei)<[BrainStateData.endtime]);
                if ~isempty(idx)
                    stateidx(timei)=find(strcmp(BrainStateData(idx).name,statenames));
                else
                    stateidx(timei)=0;
                end
            else
                stateidx(timei)=0;
            end
  
        end
        if isempty(fieldnames(v0_PSD_data))
            NEXT=1;
        else
            NEXT=length(v0_PSD_data)+1;
        end
        
        v0_PSD_data(NEXT).icV=icV;
        v0_PSD_data(NEXT).time=time;
        v0_PSD_data(NEXT).sweepnum=sweepi;
        v0_PSD_data(NEXT).PSD_max=PSD_max;
        v0_PSD_data(NEXT).PSD_mean=PSD_mean;
        v0_PSD_data(NEXT).PSD_median=PSD_median;
        v0_PSD_data(NEXT).realtime=bridgeddata(sweepi).realtime;
        v0_PSD_data(NEXT).stateidx=stateidx;
    end
end
%%
v0bins=[-90:.05:-40]/1000;
binsize=5/1000;

icV=[v0_PSD_data.icV];
time=[v0_PSD_data.time];
stateidx=[v0_PSD_data.stateidx];
PSD_max=[v0_PSD_data.PSD_max];
PSD_median=[v0_PSD_data.PSD_median];
PSD_mean=[v0_PSD_data.PSD_mean];
V0_PSD_bindata=struct;
figure(100)
clf
if isempty(statenames)
    eddig=0;
else
    eddig=length(statenames);
end
for statei=0:eddig
    if statei==0
        needed=stateidx<inf;
    else
        needed=stateidx==statei;
    end
    if ~isempty(valtozok.timeborders) & diff(valtozok.timeborders)>0
        needed=needed&time>min(valtozok.timeborders) & time<max(valtozok.timeborders);
    end
    powervals=nan(size(v0bins));
    powervalsSD=nan(size(v0bins));
    n=nan(size(v0bins));
    for bini=1:length(v0bins)
        idx=icV>=v0bins(bini)-binsize/2 & icV<v0bins(bini)+binsize/2;
        powervals(bini)=nanmean(PSD_max(idx&needed));
        powervalsSD(bini)=nanstd(PSD_max(idx&needed));
        n(bini)=sum(idx&needed);
    end
    if statei==0
        V0_PSD_bindata(statei+1).statename='All';
    else
        V0_PSD_bindata(statei+1).statename=statenames{statei};
    end
    V0_PSD_bindata(statei+1).v0bins=v0bins;
    V0_PSD_bindata(statei+1).binsize=binsize;
    V0_PSD_bindata(statei+1).powervals=powervals;
    V0_PSD_bindata(statei+1).powervalsSD=powervalsSD;
    
    subplot(length(statenames)+1,2,statei*2+1) 
    shadedErrorBar(v0bins*1000,powervals,powervalsSD)
    hold on
    plot(v0bins*1000,powervals,'k-','LineWidth',2)
    title(V0_PSD_bindata(statei+1).statename)
    xlabel('Voltage (mV)')
    ylabel('Power (V^{2})')
    xlimits=get(gca,'Xlim');
    subplot(length(statenames)+1,2,statei*2+2) 
    bar(v0bins*1000,n);
    xlim(xlimits);
    xlabel('Voltage (mV)')
    ylabel('# of points averaged')
end
    