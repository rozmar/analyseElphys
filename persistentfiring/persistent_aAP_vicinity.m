%% cutting out the vicinities of aAPs from all experiments %%% FINISH ME!!
valtozok=struct;
valtozok.radius=3; %seconds before and after the aAP
valtozok.onlysporadicAP=0;
valtozok.packets=0;
valtozok.minpacketinterval=20;%sec
valtozok.exportapstats=1;
valtozok.exportfield=1;
valtozok.exportic=1;
valtozok.exportPSD_field=1;
valtozok.sampleinterval=.001;
sortedeventdir=[dirs.eventdir,'sorted/'];

files=dir(sortedeventdir);
files([files.isdir])=[];


for xlsi=1:length(xlsdata)%:-1:
    aAPvicinity = struct;
    fieldxlsnum=find(strcmp(xlsdata(xlsi).HEKAfname,{xlsdata.HEKAfname}) & [xlsdata.field]==1);
    prevsweepnum=NaN;
    prevfieldsweepnum=NaN;
    prevPSDsweepnum=NaN;
    if  xlsdata(xlsi).field==0 & xlsdata(xlsi).juxta==0 & any(strfind(xlsdata(xlsi).anaesthesia,'awake')) & ~isempty(fieldxlsnum)
        a=dir([sortedeventdir,xlsdata(xlsi).ID,'.mat']);
        if ~isempty(a)
            load([sortedeventdir,xlsdata(xlsi).ID]);
            eventdata=persistent_sort_sporadic_persistent_aAPs(eventdata);
            if valtozok.packets==1
                %%
                aapdata=eventdata([eventdata.axonalAP]==1);
                isis=diff([-inf,[aapdata.maxtime]]);
                packetstarttimes=[aapdata((isis>valtozok.minpacketinterval)).maxtime,inf];
                apdata=struct;
                for packeti=1:length(packetstarttimes)-1
                    firstapidx=find([aapdata.maxtime]==packetstarttimes(packeti));
                    allapidx=find([aapdata.maxtime]>=packetstarttimes(packeti)&[aapdata.maxtime]<packetstarttimes(packeti+1));
                    apdata(packeti).sweepnum=aapdata(firstapidx).sweepnum;
                    apdata(packeti).maxtime=aapdata(firstapidx).maxtime;
                    apdata(packeti).maxh=aapdata(firstapidx).maxh;
                    apdata(packeti).si=aapdata(firstapidx).si;
                    apdata(packeti).baselineval=aapdata(firstapidx).baselineval;
                    apdata(packeti).apdata=aapdata(allapidx);
                end
            elseif valtozok.onlysporadicAP==1
                apdata=eventdata([eventdata.axonalAP_sporadic]==1);
            else
                apdata=eventdata([eventdata.axonalAP]==1);
            end
            if ~isempty(apdata) & ~isempty(fieldnames(apdata))
                %% intracell
                if valtozok.exportic==1
                load([dirs.bridgeddir,xlsdata(xlsi).ID],'bridgeddata');
                end
                disp([xlsdata(xlsi).ID])
                %% field
                if valtozok.exportfield==1
                    load([dirs.rawexporteddir,xlsdata(fieldxlsnum).ID],'rawdata');
                    fielddata=rawdata;
                end
                %% field PSD
                if valtozok.exportPSD_field==1
                    clear PSDdata_field
                    a=dir([dirs.PSDdir_log,xlsdata(fieldxlsnum).ID,'.mat']);
                    if ~isempty(a)
                    load([dirs.PSDdir_log,xlsdata(fieldxlsnum).ID],'PSDdata');
                    load([dirs.PSDdir_log,'stats/',xlsdata(fieldxlsnum).ID],'PSDdata_stats');
                    PSDdata_field=PSDdata;
                    clear PSDdata
                    PSDdata_field_stats=PSDdata_stats;
                    else
                        PSDdata_field=[];
                        PSDdata_field_stats=[];
                    end
                end
                for APi=1:length(apdata)
                    if isempty(fieldnames(aAPvicinity))
                        NEXT=1;
                    else
                        NEXT=length(aAPvicinity)+1;
                    end
                    sweepnum=apdata(APi).sweepnum;
                    maxtime=apdata(APi).maxtime;
                    maxh=apdata(APi).maxh;
                    
                    if valtozok.exportic==1
                        si=apdata(APi).si;
                        if prevsweepnum==sweepnum
                            y=icy;
                            downsamplerate=round(valtozok.sampleinterval/si);
                            if downsamplerate>1
                                si=valtozok.sampleinterval;
                                maxh=round(apdata(APi).maxh/downsamplerate);
                            end
                        else
                            y=bridgeddata(sweepnum).y;
                            d{1} = designfilt('lowpassiir','PassbandFrequency',2000,'StopbandFrequency',3000,'SampleRate',1/si,'DesignMethod','butter');
                            y=filtfilt(d{1},y);
                            downsamplerate=round(valtozok.sampleinterval/si);
                            if downsamplerate>1
                                y=downsample(y,downsamplerate);
                                si=valtozok.sampleinterval;
                                maxh=round(apdata(APi).maxh/downsamplerate);
                            end
                            prevsweepnum=sweepnum;
                            icy=y;
                        end
                        stepsize=round(valtozok.radius/si);
                        hossz=length(y);
                        stepback=min(stepsize,maxh-1);
                        stepforward=min(stepsize,hossz-maxh-1);
                        ic=nan(2*stepsize+1,1);
                        ic(stepsize-stepback+1:stepsize+1+stepforward)=y([-stepback:stepforward]+maxh);
                        
                        aAPvicinity(NEXT).icxlsnum=xlsi;
                        aAPvicinity(NEXT).si=si;
                        aAPvicinity(NEXT).ic=ic;
                        aAPvicinity(NEXT).timeback=stepback*si;
                        aAPvicinity(NEXT).timeforward=stepforward*si;
                    end
                    if valtozok.exportfield==1
                        fieldsweepnum=find([fielddata.realtime]==bridgeddata(sweepnum).realtime);
                        if ~isempty(fieldsweepnum)
                            si=apdata(APi).si;
                            if prevfieldsweepnum==fieldsweepnum
                                y=fieldy;
                                downsamplerate=round(valtozok.sampleinterval/si);
                                if downsamplerate>1
                                    si=valtozok.sampleinterval;
                                    maxh=round(apdata(APi).maxh/downsamplerate);
                                end
                            else
                                y=fielddata(fieldsweepnum).y;
                                d{1} = designfilt('lowpassiir','PassbandFrequency',500,'StopbandFrequency',1000,'SampleRate',1/si,'DesignMethod','butter');
                                d{2} = designfilt('highpassiir','PassbandFrequency',.4,'StopbandFrequency',.1,'SampleRate',1/si,'DesignMethod','butter');
                                y=filtfilt(d{1},y);
                                y=filtfilt(d{2},y);
                                downsamplerate=round(valtozok.sampleinterval/si);
                                if downsamplerate>1
                                    y=downsample(y,downsamplerate);
                                    si=valtozok.sampleinterval;
                                    maxh=round(apdata(APi).maxh/downsamplerate);
                                end
                                prevfieldsweepnum=fieldsweepnum;
                                fieldy=y;
                            end

                            stepsize=round(valtozok.radius/si);
                            hossz=length(y);
                            stepback=min(stepsize,maxh-1);
                            stepforward=min(stepsize,hossz-maxh-1);
                            field=nan(2*stepsize+1,1);
                            field(stepsize-stepback+1:stepsize+1+stepforward)=y([-stepback:stepforward]+maxh);
                            aAPvicinity(NEXT).field=field;
                            aAPvicinity(NEXT).time_elphys=[-stepsize:stepsize]*si;
                            aAPvicinity(NEXT).timeback=stepback*si;
                            aAPvicinity(NEXT).timeforward=stepforward*si;
                        else
                            %                             aAPvicinity(NEXT).field=[];
                            %                         aAPvicinity(NEXT).time_elphys=[-stepsize:stepsize]*si;
                            %                         aAPvicinity(NEXT).timeback=stepback*si;
                            %                         aAPvicinity(NEXT).timeforward=stepforward*si;
                        end
                    end
                    
                    if valtozok.exportPSD_field==1 & ~isempty(PSDdata_field)
                        fieldsweepnum=find([PSDdata_field.realtime]==bridgeddata(sweepnum).realtime);
                        if ~isempty(fieldsweepnum) & ~isempty(PSDdata_field(fieldsweepnum).si_powerMatrix)
                            %%
                            if prevPSDsweepnum==fieldsweepnum
                                powerMatrix=PSDy;
                                frequencyVector=frequencyVectory;
                                si_PSD=PSDdata_field(fieldsweepnum).si_powerMatrix;
                                if downsamplerate>1
                                    si_PSD=valtozok.sampleinterval;
                                end
                            else
                                powerMatrix=PSDdata_field(fieldsweepnum).powerMatrix;
                                %                             if isfield(PSDdata_field,'compress_offset')
                                %                                 powerMatrix=double(PSDdata_field(fieldsweepnum).powerMatrix)*PSDdata_field(fieldsweepnum).compress_multiplier+PSDdata_field(fieldsweepnum).compress_offset;
                                %                             end
                                frequencyVector=PSDdata_field(fieldsweepnum).frequencyVector;
                                si_PSD=PSDdata_field(fieldsweepnum).si_powerMatrix;
                                downsamplerate=round(valtozok.sampleinterval/si_PSD);
                                if downsamplerate>1
                                    powerMatrix=downsample(powerMatrix',downsamplerate)';
                                    si_PSD=valtozok.sampleinterval;
                                end
                                PSDy=powerMatrix;
                                frequencyVectory=frequencyVector;
                                prevPSDsweepnum=fieldsweepnum;
                            end
                            
                            hossz=size(powerMatrix,2);
                            stepsize=round(valtozok.radius/si_PSD);
                            maxh_field=round(apdata(APi).maxh*(apdata(APi).si/si_PSD));
                            stepback=min(stepsize,maxh_field-1);
                            stepforward=min(stepsize,hossz-maxh_field-1);
                            PSD_field=uint32(zeros(length(frequencyVector),2*stepsize+1));
                            PSD_field(:,stepsize-stepback+1:stepsize+1+stepforward)=powerMatrix(:,[-stepback:stepforward]+maxh_field);
                            aAPvicinity(NEXT).PSD_field=PSD_field;
                            aAPvicinity(NEXT).time_PSD_field=[-stepsize:stepsize]*si_PSD;
                            aAPvicinity(NEXT).frequencyVector_field=frequencyVector;
                            if isfield(PSDdata_field,'compress_offset')
                                aAPvicinity(NEXT).PSD_field_compress_multiplier=PSDdata_field(fieldsweepnum).compress_multiplier;
                                aAPvicinity(NEXT).PSD_field_compress_offset=PSDdata_field(fieldsweepnum).compress_offset;
                            else
                                aAPvicinity(NEXT).PSD_field_compress_multiplier=[];
                                aAPvicinity(NEXT).PSD_field_compress_offset=[];
                            end
%                             aAPvicinity(NEXT).timeback=stepback*si;
%                             aAPvicinity(NEXT).timeforward=stepforward*si;
                        else
                            aAPvicinity(NEXT).PSD_field=[];
                            aAPvicinity(NEXT).time_PSD_field=[];
                            aAPvicinity(NEXT).frequencyVector_field=[];
                            aAPvicinity(NEXT).PSD_field_compress_multiplier=[];
                            aAPvicinity(NEXT).PSD_field_compress_offset=[];
                        end
                    end
                    
                    
                    if valtozok.exportapstats==1
                        %%
                        types={'AP','aAP','sAP','aAP_sporadic','aAP_persistent','ep','ip'};
                        for typei=1:length(types)
                            type=types{typei};
                            eventsnow=eventdata([eventdata.maxtime]>maxtime-valtozok.radius & [eventdata.maxtime]<maxtime+valtozok.radius);
                            if strcmp(type,'AP')
                                eventsnow=eventsnow(find(strcmp({eventsnow.type},'AP')));
                            elseif strcmp(type,'ep')
                                eventsnow=eventsnow(find(strcmp({eventsnow.type},'ep')));
                            elseif strcmp(type,'ip')
                                eventsnow=eventsnow(find(strcmp({eventsnow.type},'ip')));
                            elseif strcmp(type,'aAP')
                                eventsnow=eventsnow(find([eventsnow.axonalAP]));
                            elseif strcmp(type,'sAP')
                                eventsnow=eventsnow(find([eventsnow.somaticAP]));
                            elseif strcmp(type,'aAP_sporadic')
                                eventsnow=eventsnow(find([eventsnow.axonalAP_sporadic]));
                            elseif strcmp(type,'aAP_persistent')
                                eventsnow=eventsnow(find([eventsnow.axonalAP_persistent]));
                            end
                            if isempty(eventsnow)
                                aAPvicinity(NEXT).([type,'_num'])=0;
                                aAPvicinity(NEXT).([type,'_times'])=NaN;
                                aAPvicinity(NEXT).([type,'_freq'])=0;
                                aAPvicinity(NEXT).([type,'_interval_before'])=valtozok.radius;
                                aAPvicinity(NEXT).([type,'_interval_after'])=valtozok.radius;
                            else
                                aAPvicinity(NEXT).([type,'_num'])=length(eventsnow);
                                aAPvicinity(NEXT).([type,'_times'])=[eventsnow.maxtime]-maxtime;
                                aAPvicinity(NEXT).([type,'_freq'])=length(eventsnow)/valtozok.radius/2;
                                
                                eventsbefore=eventsnow([eventsnow.maxtime]<maxtime);
                                if isempty(eventsbefore)
                                    aAPvicinity(NEXT).([type,'_interval_before'])=valtozok.radius;
                                else
                                    aAPvicinity(NEXT).([type,'_interval_before'])=abs(eventsbefore(end).maxtime-maxtime);
                                end
                                eventsafter=eventsnow([eventsnow.maxtime]>maxtime);
                                if isempty(eventsafter)
                                    aAPvicinity(NEXT).([type,'_interval_after'])=valtozok.radius;
                                else
                                    aAPvicinity(NEXT).([type,'_interval_after'])=abs(eventsafter(1).maxtime-maxtime);
                                end
                            end
                            
                        end
%                         =find
                    end
                    aAPvicinity(NEXT).baselineval=apdata(APi).baselineval;
                    aAPvicinity(NEXT).ID=xlsdata(xlsi).ID;
                     if valtozok.packets==1
                        aAPvicinity(NEXT).aAPsinpacket=apdata(APi).apdata;
                     end
                    aAPvicinity(NEXT).realtime=maxtime;
                    
                end
                if ~isempty(fieldnames(aAPvicinity))
                    if valtozok.packets==1
                        save([dirs.basedir,'Vicinity_packet/',xlsdata(xlsi).ID],'aAPvicinity','-v7.3');
                    else
                        save([dirs.basedir,'Vicinity_aAP/',xlsdata(xlsi).ID],'aAPvicinity','-v7.3');
                    end
                end
%                 disp('lol')
            end
        end
    end
end
return
%% loading the APvicinity and plotting
basedir=[dirs.basedir,'aAP_vicinity/'];
filename=uigetfile(basedir);
load([basedir,filename]);
%decompress PSD
for i=1:length(aAPvicinity)
    if isfield(aAPvicinity,'PSD_field_compress_offset') & ~isempty(aAPvicinity(i).PSD_field_compress_offset)
        aAPvicinity(i).PSD_field=double(aAPvicinity(i).PSD_field);
        aAPvicinity(i).PSD_field(aAPvicinity(i).PSD_field==0)=NaN;
        aAPvicinity(i).PSD_field=double(aAPvicinity(i).PSD_field)*aAPvicinity(i).PSD_field_compress_multiplier+aAPvicinity(i).PSD_field_compress_offset;
    end
end


figure(1)
clf
for i=1:length(aAPvicinity)
    figure(1)
    subplot(3,1,1)
    hold all
    plot(aAPvicinity(i).time_elphys,aAPvicinity(i).ic)%-aAPvicinity(i).baselineval
    
    subplot(3,1,2)
    hold all
    plot(aAPvicinity(i).time_elphys,aAPvicinity(i).field-aAPvicinity(i).field(round(length(aAPvicinity(i).field)/2)))
%     subplot(3,1,3)
% 
%     imagesc(aAPvicinity(i).time_PSD_field,aAPvicinity(i).frequencyVector_field,aAPvicinity(i).PSD_field);
%     set(gca,'YDir','normal');
%     colormap linspecer
%     ylabel('Frequency (Hz)')
%     xlabel('Time (s)')
%     
%     
%     set(gca, 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
%     ytickidxs=[];
%     yticknow=get(gca,'YTick');
%     ylimnow=get(gca,'Ylim');
%    
%     linearyticks = linspace(ylimnow(1),ylimnow(end),length(aAPvicinity(i).frequencyVector_field));
%     for ticki=1:length(yticknow)
%         [~,ytickidxs(ticki)]=min(abs(linearyticks-yticknow(ticki)));
%     end
%     yticklabelnow=aAPvicinity(i).frequencyVector_field(ytickidxs);
%     set(gca,'Yticklabel',round(yticklabelnow*100)/100)
%
%     pause
end
figure(2)
clf
subplot(3,1,1)
hold all
t=[aAPvicinity.time_elphys];
y=[aAPvicinity.ic];
plot(nanmean(t,2),nanmean(y,2))%-aAPvicinity(i).baselineval

t=[aAPvicinity.time_elphys];
y=[aAPvicinity.field];
subplot(3,1,2)
hold all
plot(nanmean(t,2),nanmean(y,2))
%%
figure(3)
clf
temp=reshape([aAPvicinity.PSD_field],[size([aAPvicinity(1).PSD_field],1),size([aAPvicinity(1).PSD_field],2),length(aAPvicinity)]);
imagesc(aAPvicinity(1).time_PSD_field,aAPvicinity(1).frequencyVector_field,nanmedian(temp,3));
set(gca,'YDir','normal');
colormap linspecer
ylabel('Frequency (Hz)')
xlabel('Time (s)')


set(gca, 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
ytickidxs=[];
yticknow=get(gca,'YTick');
ylimnow=get(gca,'Ylim');

linearyticks = linspace(ylimnow(1),ylimnow(end),length(aAPvicinity(i).frequencyVector_field));
for ticki=1:length(yticknow)
    [~,ytickidxs(ticki)]=min(abs(linearyticks-yticknow(ticki)));
end
yticklabelnow=aAPvicinity(i).frequencyVector_field(ytickidxs);
set(gca,'Yticklabel',round(yticklabelnow*100)/100)





%% 
types={'AP','aAP','sAP','aAP_sporadic','aAP_persistent','ep','ip'};
figure(3)
clf
for typei=1:length(types)
    subplot(length(types),2,typei*2-1)
    tohist=[aAPvicinity.([types{typei},'_num'])];
    hist(tohist,[0:1:1500])
    title(types{typei})
    subplot(length(types),2,typei*2)
    tohist=[aAPvicinity.([types{typei},'_times'])];
    tohist(tohist==0)=[];
    hist(tohist,[-valtozok.radius:.1:valtozok.radius])
    xlim([-valtozok.radius valtozok.radius])
end

%% 
types={'AP','aAP','sAP','aAP_sporadic','aAP_persistent','ep','ip'};
figure(4)
clf
for typei=1:length(types)
    subplot(length(types),2,typei*2-1)
    tohist=[aAPvicinity.([types{typei},'_interval_before'])];
    hist(tohist,[0:.1:10])
    title([types{typei},' before'])
    subplot(length(types),2,typei*2)
    tohist=[aAPvicinity.([types{typei},'_interval_after'])];
    tohist(tohist==0)=[];
    hist(tohist,[0:.1:10])
%     xlim([-valtozok.radius valtozok.radius])
    title([types{typei},' after'])
end
%%
figure(1)
clf
for i=1:length(aAPvicinity)
    subplot(2,1,1)
    hold all
    plot(aAPvicinity(i).time_elphys,aAPvicinity(i).ic-aAPvicinity(i).baselineval)
    
    subplot(2,1,2)
    hold all
    plot(aAPvicinity(i).time_elphys,aAPvicinity(i).field-aAPvicinity(i).field(round(length(aAPvicinity(i).field)/2)))
end
%% cutting out only the aAPs
timeback=2;
timeforward=2;
sdwindow=.5;
APwaves=struct;
IDs=unique({aAPvicinity.ID});
for IDi=1:length(IDs)
    IDnow=IDs{IDi};
    idxs=find(strcmp({aAPvicinity.ID},IDnow));
    aAPvicinity_now=aAPvicinity(idxs);
    sis=unique(round([aAPvicinity_now.si]*10^5)/10^5);
    if length(sis)>1
        disp(['different sampling intervals, resampling: ',IDnow,' - ',num2str(sis)])
        newsi=max(sis);
        for idxi=1:length(aAPvicinity_now)
            aAPvicinity_now(idxi).si = round(aAPvicinity_now(idxi).si*10^5)/10^5;
            if aAPvicinity_now(idxi).si~=newsi;
                aAPvicinity_now(idxi).ic=downsample(aAPvicinity_now(idxi).ic,round(newsi/aAPvicinity_now(idxi).si));
                aAPvicinity_now(idxi).si=newsi;
            end
        end
    end
    ic=[];
    %     icsd=[];
    for idxi=1:length(aAPvicinity_now)
        si=aAPvicinity_now(idxi).si;
        maxh=round(length(aAPvicinity_now(idxi).ic)/2);
        stepback=round(timeback/si);
        stepforward=round(timeforward/si);
        ic(:,idxi)=aAPvicinity_now(idxi).ic([-stepback:stepforward]+maxh);
        %         icsd(:,idxi)=moving(aAPvicinity_now(idxi).ic([-stepback:stepforward]+maxh),round(sdwindow/si),'std');
    end
    if isempty(fieldnames(APwaves))
        NEXT=1;
    else
        NEXT=length(APwaves)+1;
    end
    APwaves(NEXT).ic=mean(ic,2);
    %     APwaves(NEXT).icsd=mean(icsd,2);
    APwaves(NEXT).time=[-stepback:stepforward]*si;
    APwaves(NEXT).baselineval=nanmedian(APwaves(NEXT).ic(round(stepback/5)*3:round(stepback/5)*4));
    APwaves(NEXT).n=idxi;
    APwaves(NEXT).ID=IDnow;
end
%%
figure(2)
clf
hold all
for i=1:length(APwaves)
    if APwaves(i).n>1
        plot(APwaves(i).time,APwaves(i).ic-APwaves(i).baselineval);
        
    end
end