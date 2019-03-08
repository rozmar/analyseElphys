function [v0_PSD_data,valtozok] = aE_V0_vs_PSD(dirs,xlsdata,xlsnum,valtozok)
disp(['exporting PSD V0 dependence for ',xlsdata(xlsnum).ID])
V0percentilesneeded=[.1:.1:.9];
if nargin<4
    valtozok=struct;
    valtozok.movingvindowsize=3;
    valtozok.movingvindowstep=.5; %seconds for downsampling and median filtering
    valtozok.PSDonfield=true;
    valtozok.minsweeplength=valtozok.movingvindowsize/2;
end

fieldxlsnum=find(strcmp(xlsdata(xlsnum).HEKAfname,{xlsdata.HEKAfname}) & [xlsdata.field]==1);
if isempty(fieldxlsnum)
    fieldxlsnum=(NaN);
end

load([dirs.bridgeddir,xlsdata(xlsnum).ID],'bridgeddata');

if isnan(fieldxlsnum) | ~valtozok.PSDonfield
    load([dirs.PSDdir,xlsdata(xlsnum).ID]);
    PSDonfield=false;
    valtozok.PSDonfield=false;
    a=dir([dirs.PSDdir,'stats/',xlsdata(xlsnum).ID,'.mat']);
    if ~isempty(a)
        load([dirs.PSDdir,'stats/',xlsdata(xlsnum).ID,'.mat']);
    else
        PSDdata_stats=[];
    end
else
    load([dirs.PSDdir,xlsdata(fieldxlsnum).ID]);
    PSDonfield=true;
    a=dir([dirs.PSDdir,'stats/',xlsdata(fieldxlsnum).ID,'.mat']);
    if ~isempty(a)
        load([dirs.PSDdir,'stats/',xlsdata(fieldxlsnum).ID,'.mat']);
    else
        PSDdata_stats=[];
    end
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
    if si*hossz>valtozok.minsweeplength % & (isempty(valtozok.timeborders) | (sweepstart>=valtozok.timeborders(1) & sweepstart<=valtozok.timeborders(2)) | (sweepend>=valtozok.timeborders(1) & sweepend<=valtozok.timeborders(2)) | (sweepstart<=valtozok.timeborders(1) & sweepend>=valtozok.timeborders(2)))% &~isempty(fieldsweepnum)
        %%
        movingvindowstep=round(valtozok.movingvindowstep/si);
        movingvindowsize=round(valtozok.movingvindowsize/si);
        
%% single V0 median
% %         tic
% %         icV=bridgeddata(sweepi).y;
% %         time=[1:length(icV)]*si-si+bridgeddata(sweepi).realtime;
% %         %%
% %         icV=medfilt1([icV(ceil(movingvindowsize/2):-1:1),icV,icV(end:-1:end-ceil(movingvindowsize/2))],movingvindowsize);
% %         icV=icV([1:hossz]+ceil(movingvindowsize/2)+1);
% %         icV=downsample(icV,movingvindowstep);
% %         time=downsample(time,movingvindowstep);
% %         toc
% % fieldsweepnum=find([PSDdata.realtime]==bridgeddata(sweepi).realtime);
% %         si_PSD=PSDdata(fieldsweepnum).si_powerMatrix;
% %         frequencyVector=PSDdata(fieldsweepnum).frequencyVector;
% %         powerMatrix=PSDdata(fieldsweepnum).powerMatrix;
% %         neededfreqidx=frequencyVector>=min(valtozok.frequencyrange) & frequencyVector<=max(valtozok.frequencyrange);
% %         powerMatrix=powerMatrix(neededfreqidx,:);
% %         PSDtime=[1:length(PSDdata(fieldsweepnum).y)]*si_PSD+PSDdata(fieldsweepnum).realtime-si_PSD;
% %         PSD_max=zeros(size(icV));
% %         PSD_mean=zeros(size(icV));
% %         PSD_median=zeros(size(icV));
% %         stateidx=zeros(size(icV));
% %         for timei=1:length(time)
% %             needed=PSDtime>=time(timei)-valtozok.movingvindowsize/2 & PSDtime<=time(timei)+valtozok.movingvindowsize/2;
% %             datanow=powerMatrix(:,needed);
% %             PSD_max(timei)=max(datanow(:));
% %             PSD_mean(timei)=mean(datanow(:));
% %             PSD_median(timei)=median(datanow(:));
% %             
% %             if ~isempty(BrainStateData)
% %                 idx=find(time(timei)>[BrainStateData.starttime] & time(timei)<[BrainStateData.endtime]);
% %                 if ~isempty(idx)
% %                     stateidx(timei)=find(strcmp(BrainStateData(idx).name,statenames));
% %                 else
% %                     stateidx(timei)=0;
% %                 end
% %             else
% %                 stateidx(timei)=0;
% %             end
% %   
% %         end
% %         if isempty(fieldnames(v0_PSD_data))
% %             NEXT=1;
% %         else
% %             NEXT=length(v0_PSD_data)+1;
% %         end
% %         
% %         v0_PSD_data(NEXT).icV=icV;
% %         v0_PSD_data(NEXT).time=time;
% %         v0_PSD_data(NEXT).sweepnum=sweepi;
% %         v0_PSD_data(NEXT).PSD_max=PSD_max;
% %         v0_PSD_data(NEXT).PSD_mean=PSD_mean;
% %         v0_PSD_data(NEXT).PSD_median=PSD_median;
% %         v0_PSD_data(NEXT).realtime=bridgeddata(sweepi).realtime;
% %         v0_PSD_data(NEXT).stateidx=stateidx;
        %%

        icV=bridgeddata(sweepi).y;
        time_orig=[1:length(icV)]*si-si+bridgeddata(sweepi).realtime;

        icVs=zeros(length(V0percentilesneeded),ceil(hossz/movingvindowstep));
        time=[];
        for stepi=1:size(icVs,2)
            sortedy=sort(icV((stepi-1)*movingvindowstep+1:min(length(icV),(stepi-1)*movingvindowstep+movingvindowsize)));
            icVs(:,stepi)=sortedy(ceil(V0percentilesneeded*length(sortedy)));
            time(stepi)=time_orig(round(mean([(stepi-1)*movingvindowstep+1,min(length(icV),(stepi-1)*movingvindowstep+movingvindowstep)])));
        end

        fieldsweepnum=find([PSDdata.realtime]==bridgeddata(sweepi).realtime);
        si_PSD=PSDdata(fieldsweepnum).si_powerMatrix;
        frequencyVector=PSDdata(fieldsweepnum).frequencyVector;
        if isfield(PSDdata,'compress_offset')
%             powerMatrix=double(PSDdata(fieldsweepnum).powerMatrix);
            powerMatrix=double(PSDdata(fieldsweepnum).powerMatrix)*PSDdata(fieldsweepnum).compress_multiplier+PSDdata(fieldsweepnum).compress_offset;
        else
            powerMatrix=PSDdata(fieldsweepnum).powerMatrix;
        end
        
%         neededfreqidx=frequencyVector>=min(valtozok.frequencyrange) & frequencyVector<=max(valtozok.frequencyrange);
%         powerMatrix=powerMatrix(neededfreqidx,:);
        PSDtime=[1:length(PSDdata(fieldsweepnum).y)]*si_PSD+PSDdata(fieldsweepnum).realtime-si_PSD;
        PSD_max=zeros(length(frequencyVector),size(icVs,2));
        PSD_mean=zeros(length(frequencyVector),size(icVs,2));
        PSD_median=zeros(length(frequencyVector),size(icVs,2));
        stateidx=zeros(1,size(icVs,2));
        %%
        for timei=1:length(time)
            needed=PSDtime>=time(timei)-valtozok.movingvindowsize/2 & PSDtime<=time(timei)+valtozok.movingvindowsize/2;
            datanow=powerMatrix(:,needed);
            PSD_max(:,timei)=max(datanow,[],2);
            PSD_mean(:,timei)=mean(datanow,2);
            PSD_median(:,timei)=median(datanow,2);
            
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
        
        v0_PSD_data(NEXT).icVs=icVs;
        v0_PSD_data(NEXT).time=time;
        v0_PSD_data(NEXT).sweepnum=sweepi;
        v0_PSD_data(NEXT).PSD_max=PSD_max;
        v0_PSD_data(NEXT).PSD_mean=PSD_mean;
        v0_PSD_data(NEXT).PSD_median=PSD_median;
        v0_PSD_data(NEXT).frequencyVector=frequencyVector;
        v0_PSD_data(NEXT).v0_percentiles=V0percentilesneeded;
        v0_PSD_data(NEXT).realtime=bridgeddata(sweepi).realtime;
        v0_PSD_data(NEXT).stateidx=stateidx;
        v0_PSD_data(NEXT).statenames=statenames;
    end
end

% persistent_V0_vs_PSD_GUI(v0_PSD_data,dirs,xlsdata,xlsnum,valtozok)

% %%
% percentileidx=1;
% 
% v0bins=[-90:.01:-40]/1000;
% binsize=1/1000;
% 
% icVs=[v0_PSD_data.icVs];
% icV=icVs(percentileidx,:);
% time=[v0_PSD_data.time];
% stateidx=[v0_PSD_data.stateidx];
% PSD_max=[v0_PSD_data.PSD_max];
% PSD_median=[v0_PSD_data.PSD_median];
% PSD_mean=[v0_PSD_data.PSD_mean];
% PSD_to_use=zscore(PSD_median,0,2);
% % PSD_to_use=PSD_max;
% V0_PSD_bindata=struct;
% frequencyVector=v0_PSD_data(1).frequencyVector;
% 
% %%
% figure(100)
% clf
% if isempty(statenames)
%     eddig=0;
% else
%     eddig=length(statenames);
% end
% for statei=0:eddig
%     if statei==0
%         needed=stateidx<inf;
%     else
%         needed=stateidx==statei;
%     end
%     if ~isempty(valtozok.timeborders) & diff(valtozok.timeborders)>0
%         needed=needed&time>min(valtozok.timeborders) & time<max(valtozok.timeborders);
%     end
%     powervals=nan(length(frequencyVector),length(v0bins));
%     powervalsSD=nan(length(frequencyVector),length(v0bins));
%     n=nan(size(v0bins));
%     for bini=1:length(v0bins)
%         idx=icV>=v0bins(bini)-binsize/2 & icV<v0bins(bini)+binsize/2;
%         powervals(:,bini)=nanmean(PSD_to_use(:,idx&needed),2);
%         powervalsSD(:,bini)=nanstd(PSD_to_use(:,idx&needed),1,2);
%         n(bini)=sum(idx&needed);
%     end
%     if statei==0
%         V0_PSD_bindata(statei+1).statename='All';
%     else
%         V0_PSD_bindata(statei+1).statename=statenames{statei};
%     end
%     V0_PSD_bindata(statei+1).v0bins=v0bins;
%     V0_PSD_bindata(statei+1).binsize=binsize;
%     V0_PSD_bindata(statei+1).powervals=powervals;
%     V0_PSD_bindata(statei+1).powervalsSD=powervalsSD;
%     
%     subplot(length(statenames)+1,3,statei*3+1) 
% %     shadedErrorBar(v0bins*1000,powervals,powervalsSD)
%     imagesc(v0bins*1000,frequencyVector,powervals)
%     set(gca,'YDir','normal');
%     colormap linspecer
%     
%     xlabel('Voltage (mV)')
%     ylabel('Frequency (Hz)')
% %     ylabel('Power (V^{2})')
%     xlimits=get(gca,'Xlim');
%     title('mean')
%     subplot(length(statenames)+1,3,statei*3+2) 
% %     shadedErrorBar(v0bins*1000,powervals,powervalsSD)
%     imagesc(v0bins*1000,frequencyVector,powervalsSD)
%     set(gca,'YDir','normal');
%     colormap linspecer
%     
%     xlabel('Voltage (mV)')
%     ylabel('Frequency (Hz)')
% %     ylabel('Power (V^{2})')
%     xlimits=get(gca,'Xlim');
%     title('sd')
%     subplot(length(statenames)+1,3,statei*3+3) 
%     bar(v0bins*1000,n);
%     xlim(xlimits);
%     xlabel('Voltage (mV)')
%     ylabel('# of points averaged')
%     title(V0_PSD_bindata(statei+1).statename)
% end
%     