%% cutting out the vicinities of aAPs from all experiments %%% FINISH ME!!
valtozok=struct;
valtozok.radius=10; %seconds before and after the aAP
valtozok.onlysporadicAP=1;
valtozok.exportapstats=1;
valtozok.exportfield=0;
valtozok.exportic=0;

sortedeventdir=[dirs.eventdir,'sorted/'];

files=dir(sortedeventdir);
files([files.isdir])=[];
aAPvicinity = struct;

for xlsi=1:length(xlsdata)
    fieldxlsnum=find(strcmp(xlsdata(xlsi).HEKAfname,{xlsdata.HEKAfname}) & [xlsdata.field]==1);
    if  xlsdata(xlsi).field==0 & xlsdata(xlsi).juxta==0 & any(strfind(xlsdata(xlsi).anaesthesia,'awake')) & ~isempty(fieldxlsnum)
        a=dir([sortedeventdir,xlsdata(xlsi).ID,'.mat']);
        if ~isempty(a)
            load([sortedeventdir,xlsdata(xlsi).ID]);
            eventdata=persistent_sort_sporadic_persistent_aAPs(eventdata);
            if valtozok.onlysporadicAP==1
                apdata=eventdata([eventdata.axonalAP_sporadic]==1);
            else
                apdata=eventdata([eventdata.axonalAP]==1);
            end
            if ~isempty(apdata)>0
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
                %%
                for APi=1:length(apdata)
                    %%
                    if isempty(fieldnames(aAPvicinity))
                        NEXT=1;
                    else
                        NEXT=length(aAPvicinity)+1;
                    end
                    sweepnum=apdata(APi).sweepnum;
                    maxtime=apdata(APi).maxtime;
                    maxh=apdata(APi).maxh;
                    si=apdata(APi).si;
                    
                    stepsize=round(valtozok.radius/si);
                    
                    if valtozok.exportic==1
                        
                        y=bridgeddata(sweepnum).y;
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
                        y=fielddata(fieldsweepnum).y;
                        hossz=length(y);
                        stepback=min(stepsize,maxh-1);
                        stepforward=min(stepsize,hossz-maxh-1);
                        field=nan(2*stepsize+1,1);
                        field(stepsize-stepback+1:stepsize+1+stepforward)=y([-stepback:stepforward]+maxh);
                        aAPvicinity(NEXT).field=field;
                        aAPvicinity(NEXT).time_elphys=[-stepsize:stepsize]*si;
                        aAPvicinity(NEXT).timeback=stepback*si;
                        aAPvicinity(NEXT).timeforward=stepforward*si;
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
                    
                    
                end
                disp('lol')
            end
        end
    end
end
return
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