%%
close all
clear all
projectnum=2;

%%
% alapadatok


if projectnum==1;
    overwrite=0;
    locations=marcicucca_locations;
    dirs.basedir=[locations.tgtardir,'ANALYSISdata/marci/Human_rosehip/CB1elphys/'];
    dirs.rawexporteddir=[dirs.basedir,'Exported_raw/'];
    dirs.bridgeddir=[dirs.basedir,'Bridged_stim/'];
    dirs.eventdir=[dirs.basedir,'Events/'];
    dirs.onlyAPeventdir=[dirs.basedir,'Events_onlyAP/'];
    dirs.grpupedeventdir=[dirs.basedir,'Events_grouped/'];
    dirs.stimepochdir=[dirs.basedir,'Stimepochs/'];
    dirs.figuresdir=[dirs.basedir,'Figures/'];
    amplifier='HEKA';
    xlsdata=aE_readxls([dirs.basedir,'cb1elphys.xls']);
elseif projectnum==2;
    overwrite=0;
    locations=marcicucca_locations;
    dirs.basedir=[locations.tgtardir,'ANALYSISdata/marci/_persistent/_InVivo/'];
    dirs.rawdir=[locations.tgtardir,'AXONdata/'];
    dirs.bridgeddir=[dirs.basedir,'Bridged_stim/'];
    dirs.eventdir=[dirs.basedir,'Events/'];
    % dirs.onlyAPeventdir=[dirs.basedir,'Events_onlyAP/'];
    % dirs.grpupedeventdir=[dirs.basedir,'Events_grouped/'];
    % dirs.stimepochdir=[dirs.basedir,'Stimepochs/'];
    % dirs.figuresdir=[dirs.basedir,'Figures/'];
    xlsdata=aE_readxls([dirs.basedir,'invivodata.xls']);
    amplifier='AXON';
end

if strcmp(amplifier,'AXON')
    aE_exportAXONdata(dirs,xlsdata,overwrite)
else
    %% export Raw data from HEKA
    analysElphys_exportrawHEKAdata(dirs,xlsdata,overwrite)
    %% generate PGF data, bridge balancing
    aE_generatePGF_bridge_HEKA(dirs,xlsdata,overwrite)
end

%% finding events
valtozok.overwrite=1;
valtozok.plotit=0;
valtozok.threshholdaveragetime=30;%s
% valtozok.mindvpdt=1;
valtozok.minampl=.0001;
valtozok.apampl=.005;
% valtozok.maxfilteredhw=.01;
valtozok.maxdecayriseratio=5;
valtozok.maxrisetime=.01;
valtozok.minaphw=.0001;
valtozok.maxapwidth=.005;
valtozok.filtermovingtime=.0025;
valtozok.diffmovingt=.0005;
valtozok.steptime=.0005; %s
valtozok.eventminsdval=3;
valtozok.apthreshval=10;
aE_findevents(valtozok,dirs)

% return
%%  inspecting spikes
% [Selection,ok] = listdlg('ListString',{xlsdata.ID},'ListSize',[300 600]); % az XLS file alapján kiválasztjuk, hogy melyik file összes mérésén szeretnénk végigmenni
Selection=1;
APdata=struct;
timebefore=.002;
timeafter=.004;
movingaverage=4;
for xlsnum=1:length(Selection) %going throught potential presynaptic cells
    preevents=load([dirs.eventdir,xlsdata(prenum).ID]);
    pretraces=load([dirs.bridgeddir,xlsdata(prenum).ID]);
    preevents.eventdata=preevents.eventdata(strcmp({preevents.eventdata.type},'AP'));%&[preevents.eventdata.amplitude]>0.03
    
    apdiffs=diff([preevents.eventdata.maxtime]);
    preevents.eventdata(find(apdiffs==0)+1)=[];
    apdiffs=diff([preevents.eventdata.maxtime]);
    prediffs=[inf,apdiffs];
    postdiffs=[apdiffs,inf];
    postpostdiffs=[apdiffs(2:end),inf,inf]; %time after the second spike in the paired pulse protocol should be more than this value
%     neededapidx=find(prediffs>1 & postdiffs<.065 & postdiffs>.045 & postpostdiffs>.1); %criteria to select presynaptic APs
    neededapidx=find([preevents.eventdata.amplitude]>0.0);
    neededapmaxtimes=[preevents.eventdata(neededapidx).maxtime];
    neededapmaxhs=[preevents.eventdata(neededapidx).maxh];
    neededapsweepnums=[preevents.eventdata(neededapidx).sweepnum];
    for APnumireal=1:length(neededapidx)
        if isempty(fieldnames(APdata))
            APnumi=1;
        else
            APnumi=length(APdata)+1;
        end
        sweepnum=preevents.eventdata(neededapidx(APnumireal)).sweepnum;
        maxh=preevents.eventdata(neededapidx(APnumireal)).maxh;
        si=pretraces.bridgeddata(sweepnum).si;
        stepback=round(timebefore/si);
        stepforward=round(timeafter/si);
        if maxh-stepback*2>0 & maxh+stepforward*2<length(pretraces.bridgeddata(sweepnum).y)
            apwave=moving(pretraces.bridgeddata(sweepnum).y(maxh-stepback:maxh),movingaverage);
            dapwave=diff(apwave)/si;
            [~,dvmaxh]=max(dapwave);
            threshh=find(dapwave(1:dvmaxh)<10,1,'last');
            dddapwave=(diff(apwave(threshh:end),2));
            [~,maxdddhely]=max(dddapwave);
            
            difi=abs(threshh-stepback);
            maxh=maxh-difi;
            if ~isempty(maxh)
            APdata(APnumi).y=pretraces.bridgeddata(sweepnum).y(maxh-stepback:maxh+stepforward)';
            APdata(APnumi).x=[-stepback*si:si:stepforward*si]';
            APdata(APnumi).sweepnum=sweepnum;
            APdata(APnumi).starth=maxh-stepback;
            APdata(APnumi).endh=maxh+stepforward;
            
            APdata(APnumi).si=si;
            APdata(APnumi).dVperdt=diff(moving(APdata(APnumi).y,movingaverage))/si;
            APdata(APnumi).dVperdtT=APdata(APnumi).x(1:end-1)+si/2;
            APdata(APnumi).dVperdtV=mean([APdata(APnumi).y(1:end-1),APdata(APnumi).y(2:end)],2);
            
            APdata(APnumi).ddVperdt=diff(APdata(APnumi).dVperdt)/si;
            APdata(APnumi).ddVperdtT=APdata(APnumi).dVperdtT(1:end-1)+si/2;
            APdata(APnumi).ddVperdtV=mean([ APdata(APnumi).dVperdtV(1:end-1), APdata(APnumi).dVperdtV(2:end)],2);
            
            APdata(APnumi).dddVperdt=diff(APdata(APnumi).ddVperdt)/si;
            
            APdata(APnumi).dddVperdtT=APdata(APnumi).ddVperdtT(1:end-1)+si/2;
            APdata(APnumi).dddVperdtV=mean([ APdata(APnumi).ddVperdtV(1:end-1), APdata(APnumi).ddVperdtV(2:end)],2);
            
            if min(dddapwave(1:maxdddhely))<0
                APdata(APnumi).persistent=1;
%                 figure(1)
%                 clf
%                 plot(dddapwave)
% %                 title()
%                 pause
            else
                APdata(APnumi).persistent=0;
            end
            end
            
        end
    end 
end
%
APdata=APdata([APdata.si]==min([APdata.si]));
APdataold=APdata;

APdata=APdataold([APdataold.persistent]==0);
%
% close all
figure(2)
clf
h1=subplot(2,4,1)
plot([APdata.x],[APdata.y])
hold on
plot(mean([APdata.x],2),mean([APdata.y],2),'k-','LineWidth',3)
axis tight
h2=subplot(2,4,2)
plot([APdata.dVperdtT],[APdata.dVperdt])
hold on
plot(mean([APdata.dVperdtT],2),mean([APdata.dVperdt],2),'k-','LineWidth',3)
axis tight
h3=subplot(2,4,3)
plot([APdata.ddVperdtT],[APdata.ddVperdt])
hold on
plot(mean([APdata.ddVperdtT],2),mean([APdata.ddVperdt],2),'k-','LineWidth',3)
axis tight
h4=subplot(2,4,4)
plot([APdata.dVperdtV], [APdata.dVperdt])
hold on
plot(mean([APdata.dVperdtV],2),mean([APdata.dVperdt],2),'k-','LineWidth',3)
axis tight

APdata=APdataold([APdataold.persistent]==1);
h5=subplot(2,4,5)
plot([APdata.x],[APdata.y])
hold on
plot(mean([APdata.x],2),mean([APdata.y],2),'k-','LineWidth',3)
axis tight
h6=subplot(2,4,6)
plot([APdata.dVperdtT],[APdata.dVperdt])
hold on
plot(mean([APdata.dVperdtT],2),mean([APdata.dVperdt],2),'k-','LineWidth',3)
axis tight
h7=subplot(2,4,7)
plot([APdata.ddVperdtT],[APdata.ddVperdt])
hold on
plot(mean([APdata.ddVperdtT],2),mean([APdata.ddVperdt],2),'k-','LineWidth',3)
axis tight
h8=subplot(2,4,8)
plot([APdata.dVperdtV], [APdata.dVperdt])
hold on
plot(mean([APdata.dVperdtV],2),mean([APdata.dVperdt],2),'k-','LineWidth',3)
axis tight
linkaxes([h1,h5],'xy')
linkaxes([h2,h6],'xy')
linkaxes([h3,h7],'xy')
linkaxes([h4,h8],'xy')
%
figure(3)
clf
hold on
for sweepnum=1:length(pretraces.bridgeddata)
    idxesnow=find([APdata.sweepnum]==sweepnum);
    if ~isempty(idxesnow)
    si=pretraces.bridgeddata(sweepnum).si;
    hossz=length(pretraces.bridgeddata(sweepnum).y);
    time=[0:si:((hossz-1)*si)]+pretraces.bridgeddata(sweepnum).realtime;
    plot(time,pretraces.bridgeddata(sweepnum).y,'k-')
    
    for i=1:length(idxesnow)
        hs=[APdata(idxesnow(i)).starth:APdata(idxesnow(i)).endh];
        plot(time(hs),pretraces.bridgeddata(sweepnum).y(hs),'r-','LineWidth',2)
    end
    end
end
return
%% check electrotonic and chemical connectivity 

dpi=600;
xcm=20;
ycm=14;

betumeret=8;
axesvastagsag=2;


xinch=xcm/2.54;
yinch=ycm/2.54;
xsize=dpi*xinch;
ysize=dpi*yinch;

clear valtozok
valtozok.gj_baselinelength=.010;
valtozok.gj_baselinelengthend=.08;
valtozok.gj_minlinelength=.05;
valtozok.gj_mincurrampl=-10*10^-12;
valtozok.baselinelength=0.025; %s
valtozok.psplength=.15; %s
valtozok.filterorder=3;
valtozok.cutofffreq=1500;
valtozok.drugwashintime=120;
valtozok.maxy0baselinedifference=.0005;
valtozok.discardpostsweepswithap=1;

[Selection,ok] = listdlg('ListString',{xlsdata.ID},'ListSize',[300 600]); % az XLS file alapján kiválasztjuk, hogy melyik file összes mérésén szeretnénk végigmenni

for xlsnum=1:length(Selection) %going throught potential presynaptic cells
    prenum=Selection(xlsnum);
    disp(['potential presynaptic cell: ',xlsdata(prenum).ID]);
    potpostidx=strcmp(xlsdata(prenum).HEKAfname,{xlsdata.HEKAfname}); %selecting the cells that may have been recorded simultaneously
    potpostidx(prenum)=0;
    dirs.figuresdirnow=[dirs.figuresdir,xlsdata(prenum).ID,'/'];
    a=dir(dirs.figuresdirnow);
    if isempty(a)
        mkdir(dirs.figuresdirnow);
    end
    if sum(potpostidx)>0
        preevents=load([dirs.eventdir,xlsdata(prenum).ID]);
        pretraces=load([dirs.bridgeddir,xlsdata(prenum).ID]);
        preevents.eventdata=preevents.eventdata(strcmp({preevents.eventdata.type},'AP'));
        apdiffs=diff([preevents.eventdata.maxtime]);
        preevents.eventdata(find(apdiffs==0)+1)=[];
        apdiffs=diff([preevents.eventdata.maxtime]);
        prediffs=[inf,apdiffs];
        postdiffs=[apdiffs,inf];
        postpostdiffs=[apdiffs(2:end),inf,inf]; %time after the second spike in the paired pulse protocol should be more than this value
        neededapidx=find(prediffs>1 & postdiffs<.065 & postdiffs>.045 & postpostdiffs>.1); %criteria to select presynaptic APs
        neededapmaxtimes=[preevents.eventdata(neededapidx).maxtime];
        neededapmaxhs=[preevents.eventdata(neededapidx).maxh];
        neededapsweepnums=[preevents.eventdata(neededapidx).sweepnum];
        findpostidxes=find(potpostidx);
        %% finding hyperpolarizing steps for GJ testing
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
%             figure(1)
%             clf
%             plot(pretraces.stimdata(sweepnum).y)
%             title([num2str(tdiffs*1000),' - ',num2str(currdiffs*10^12) ])
%             pause
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
        
        %%
        for potpostnum=1:length(findpostidxes) %going throught potential postsynaptic cells
            tracedata=struct;
            tracedataGJ=struct;
            posttraces=load([dirs.bridgeddir,xlsdata(findpostidxes(potpostnum)).ID]);
            
            if valtozok.discardpostsweepswithap==1
                postevents=load([dirs.eventdir,xlsdata(findpostidxes(potpostnum)).ID]);
                postevents.eventdata=postevents.eventdata(strcmp({postevents.eventdata.type},'AP'));
                apdiffs=diff([postevents.eventdata.maxtime]);
                postevents.eventdata(find(apdiffs==0)+1)=[];
            else
                postevents=preevents;
                postevents.eventdata=postevents.eventdata(end);
            end
            
            
            if ~isempty(fieldnames(prehypsweeps))
            for prehypsweepnum=1:length(prehypsweeps)
                postsweepnum=find(prehypsweeps(prehypsweepnum).realtime==[posttraces.bridgeddata.realtime] | prehypsweeps(prehypsweepnum).timertime==[posttraces.bridgeddata.timertime]) ;
%                  if ~isempty(postsweepnum) 
%                 figure(22)
%                 clf
%                 subplot(2,2,1)
%                 plot(pretraces.bridgeddata(prehypsweeps(prehypsweepnum).sweepnum).y)
%                 title('pre')
%                 subplot(2,2,3)
%                 plot(pretraces.stimdata(prehypsweeps(prehypsweepnum).sweepnum).y)
%                 subplot(2,2,2)
%                 plot(posttraces.bridgeddata(postsweepnum).y)
%                 title('post')
%                 subplot(2,2,4)
%                 plot(posttraces.stimdata(postsweepnum).y)
%                 pause
%                  end
                if ~isempty(postsweepnum)  & ~any([postevents.eventdata.maxtime]>prehypsweeps(prehypsweepnum).realtime+prehypsweeps(prehypsweepnum).starth*prehypsweeps(prehypsweepnum).si-valtozok.gj_baselinelength & [postevents.eventdata.maxtime]<prehypsweeps(prehypsweepnum).realtime+prehypsweeps(prehypsweepnum).endh*prehypsweeps(prehypsweepnum).si+valtozok.gj_baselinelength)% 
                    % csak akkor érdekel minket, hogyha nincs beinjektált
                    % áram a másik sejtben
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
%                         figure(1234)
%                         clf
%                         subplot(4,1,1)
%                         plot(tracedataGJ(NEXT).pre_y)
%                         subplot(4,1,2)
%                         plot(pretraces.stimdata(prehypsweeps(prehypsweepnum).sweepnum).y)
%                         
%                         subplot(4,1,3)
%                         plot(tracedataGJ(NEXT).post_y)
%                         subplot(4,1,4)
%                         plot(poststim)
%                         title(num2str(tracedataGJ(NEXT).post_realtime))
% %                         return
% %                         pause
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
%                 return
                %%
                for sweepnum=1:length(tracedataGJ)
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
                if length(neededtraces)>1
                    time=tracedataGJ(neededtraces(1)).time;
                    %plot GJ coupling
                    figure(1233)
                    clf
                    hold on
                    zeroed_pre_y=bsxfun(@(x,y) x-y, [tracedataGJ(neededtraces).pre_y], [tracedataGJ(neededtraces).pre_y0]);
                    zeroed_post_y=bsxfun(@(x,y) x-y, [tracedataGJ(neededtraces).post_y], [tracedataGJ(neededtraces).post_y0]);
                    plot(time*1000,zeroed_pre_y'*1000,'k-','Color',[.8 .8 .8])
                    plot(time*1000,nanmean(zeroed_pre_y')*1000,'k-','LineWidth',2)
                    axis tight
                    ylimits=ylim;
                    maxval=max(nanmean(zeroed_pre_y'));
                    minval=min(nanmean(zeroed_pre_y'));
                    dval=(maxval-minval)*1000;
                    ylimits(1)=max(ylimits(1),-dval*3);
                    ylimits(2)=min(ylimits(2),dval*3);
                    ylim(ylimits)
%                     return
                    xlabel('time (ms)')
                    ylabel('pre Vm (mV)')
                    set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Position',[1/xcm 1/ycm 1-2/xcm 1-2/ycm])%,'Xtick',[-valtozok.baselinelength:valtozok.baselinelength:valtozok.psplength]*1000
                    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 xsize/dpi ysize/dpi])
                    print(gcf,[dirs.figuresdirnow,xlsdata(prenum).ID,'-to-',xlsdata(findpostidxes(potpostnum)).ID,'-GJ-pre.jpg'],'-djpeg',['-r',num2str(dpi)])
                    figure(1233)
                    clf
                    hold on
                    plot(time*1000,zeroed_post_y'*1000,'k-','Color',[.8 .8 .8])
                    plot(time*1000,nanmean(zeroed_post_y')*1000,'k-','LineWidth',2)
                    axis tight
                    ylimits=ylim;
                    maxval=max(nanmean(zeroed_post_y'));
                    minval=min(nanmean(zeroed_post_y'));
                    dval=(maxval-minval)*1000;
                    ylimits(1)=max(ylimits(1),-dval*3);
                    ylimits(2)=min(ylimits(2),dval*3);
                    ylim(ylimits)
%                     return
                    xlabel('time (ms)')
                    ylabel('pre Vm (mV)')
                    set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Position',[1/xcm 1/ycm 1-2/xcm 1-2/ycm])%,'Xtick',[-valtozok.baselinelength:valtozok.baselinelength:valtozok.psplength]*1000
                    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 xsize/dpi ysize/dpi])
                    print(gcf,[dirs.figuresdirnow,xlsdata(prenum).ID,'-to-',xlsdata(findpostidxes(potpostnum)).ID,'-GJ-post.jpg'],'-djpeg',['-r',num2str(dpi)])
                end
                
                
                
            end
            for preapnum=1:length(neededapmaxtimes)
                postsweepidx=find(pretraces.bridgeddata(neededapsweepnums(preapnum)).realtime==[posttraces.bridgeddata.realtime] | pretraces.bridgeddata(neededapsweepnums(preapnum)).timertime==[posttraces.bridgeddata.timertime]);                
                if ~isempty(postsweepidx) & ~any([postevents.eventdata.maxtime]>neededapmaxtimes(preapnum)-valtozok.baselinelength & [postevents.eventdata.maxtime]<neededapmaxtimes(preapnum)+valtozok.psplength)% 
                    if isempty(fieldnames(tracedata))
                        NEXT=1;
                    else
                        NEXT=length(tracedata)+1;
                    end
                    fieldek=fieldnames(pretraces.bridgeddata);
                    for finum=1:length(fieldek)
                        tracedata(NEXT).(['pre_',fieldek{finum}])=pretraces.bridgeddata(neededapsweepnums(preapnum)).(fieldek{finum});
                        tracedata(NEXT).(['post_',fieldek{finum}])=posttraces.bridgeddata(postsweepidx).(fieldek{finum});
                        tracedata(NEXT).apmaxh=neededapmaxhs(preapnum);
                        
                    end
                end
%                 figure(1)
%                 clf
%                         subplot(2,1,1)
%                         plot(pretraces.bridgeddata(neededapsweepnums(preapnum)).y)
%                         subplot(2,1,2)
%                         plot(pretraces.stimdata(neededapsweepnums(preapnum)).y)
%                         title(num2str([preevents.eventdata(neededapidx(preapnum)).injectedcurrentatpeak,preevents.eventdata(neededapidx(preapnum)).injectedrelativecurrentatonset,preevents.eventdata(neededapidx(preapnum)).stimulated]))
%                         pause
            end
%             return
            if ~isempty(fieldnames(tracedata))
            [si,ia,ic]=unique(round([tracedata.post_si]*10^6));
            si=si/10^6;
            if length(si)>1
                disp(['si gebasz']) % if not all sweeps are recorded with the same sampling freqency, there might be some problem
                si=max(si);
%                 return
                %%
                for sweepnum=1:length(tracedata)
                    if round((tracedata(sweepnum).pre_si*10^6))/10^6<si
                        tracedata(sweepnum).apmaxh=ceil(tracedata(sweepnum).apmaxh*round((tracedata(sweepnum).pre_si*10^6))/(si*10^6));
                        tracedata(sweepnum).pre_y=resample(tracedata(sweepnum).pre_y,round((tracedata(sweepnum).pre_si*10^6)),si*10^6);
                        tracedata(sweepnum).post_y=resample(tracedata(sweepnum).post_y,round((tracedata(sweepnum).pre_si*10^6)),si*10^6);
                        tracedata(sweepnum).post_si=si;
                        tracedata(sweepnum).pre_si=si;
%                         figure(1)
%                         clf
%                         plot(tracedata(sweepnum).pre_y)
%                         hold on
%                         plot(tracedata(sweepnum).apmaxh,tracedata(sweepnum).pre_y(tracedata(sweepnum).apmaxh),'ro')
%                         pause
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
                    
%                     tracedata(sweepnum).pre_y=tracedata(sweepnum).pre_y(tracedata(sweepnum).apmaxh-stepback:tracedata(sweepnum).apmaxh+stepforward);
%                     tracedata(sweepnum).post_y=tracedata(sweepnum).post_y(tracedata(sweepnum).apmaxh-stepback:tracedata(sweepnum).apmaxh+stepforward);

                    tracedata(sweepnum).post_y0=nanmean(tracedata(sweepnum).post_y(1:stepback));
                    tracedata(sweepnum).pre_y0=nanmean(tracedata(sweepnum).pre_y(1:stepback));
                    tracedata(sweepnum).post_y0sd=nanstd(tracedata(sweepnum).post_y(1:stepback));
                
                   
%                     figure(333)

%                     sweeptodel=[sweeptodel,sweepnum];
            end
%             tracedata(sweeptodel)=[];
            % plotting paired pulses
            % if we go past 0 AM during the recording, there may be some
            % problem
%             if pretraces.bridgeddata([pretraces.bridgeddata.realtime]<pretraces.bridgeddata.realtime(1)).realtime=pretraces.bridgeddata([pretraces.bridgeddata.realtime]<pretraces.bridgeddata.realtime(1)).realtime+24*3600;
%                if [pretraces.bridgeddata.realtime] >24*3600
%                    pause
%                end
            if xlsdata(prenum).drugnum>0
                ctrlidxes=find([tracedata.pre_realtime]<min([xlsdata(prenum).drugdata.DrugWashinTime]));
            else
                ctrlidxes=find([tracedata.pre_realtime]);
            end
            y0s=[tracedata.post_y0];
            prey0s=[tracedata.pre_y0];
            y0stds=[tracedata.post_y0sd];
            zeroed_post_y=bsxfun(@(x,y) x-y, [tracedata.post_y], y0s);
                            % removing noisy sweeps

            voltmar=0;
                while voltmar==0 | any(median(ctrldiffs)+3*std(ctrldiffs)<ctrldiffs) | any(median(ctrldiffs)-3*std(ctrldiffs)>ctrldiffs)
                    ctrldiffs=bsxfun(@(x,y) x-y, [zeroed_post_y(:,ctrlidxes)], mean(zeroed_post_y(:,ctrlidxes),2));
                    ctrldiffs=max(abs(ctrldiffs));
                    if any(median(ctrldiffs)+3*std(ctrldiffs)<ctrldiffs)
                        [~,idxtodel]=max(ctrldiffs);
                        ctrlidxes(idxtodel)=[];
                    end
                    if any(median(ctrldiffs)-3*std(ctrldiffs)>ctrldiffs)
                        [~,idxtodel]=min(ctrldiffs);
                        ctrlidxes(idxtodel)=[];
                    end
                    voltmar=1;
                end
              % deleting outlying y0 values
                while any(median(y0s(ctrlidxes))+3*std(y0s(ctrlidxes))<y0s(ctrlidxes)) | any(median(y0s(ctrlidxes))-3*std(y0s(ctrlidxes))>y0s(ctrlidxes))
                    if any(median(y0s(ctrlidxes))+3*std(y0s(ctrlidxes))<y0s(ctrlidxes))
                        [~,idxtodel]=max(y0s(ctrlidxes));
                        ctrlidxes(idxtodel)=[];
                    end
                    if any(median(y0s(ctrlidxes))-3*std(y0s(ctrlidxes))>y0s(ctrlidxes))
                        [~,idxtodel]=min(y0s(ctrlidxes));
                        ctrlidxes(idxtodel)=[];
                    end
                end  
                 figure(1)
                clf
                
%                 subplot(2,1,1)
                hold on
%                 plot(time*1000,(bsxfun(@(x,y) x-y, [tracedata(ctrlidxes).pre_y], prey0s(ctrlidxes))')*1000,'k-','Color',[.8 .8 .8])
%                 plot(time*1000,(bsxfun(@(x,y) x-y, [tracedata(drugidxes).pre_y], prey0s(drugidxes))')*1000,'r-','Color',[.9 .4 .4])
%                 plot(time*1000,mean(bsxfun(@(x,y) x-y, [tracedata(ctrlidxes).pre_y], prey0s(ctrlidxes))')*1000,'k-','LineWidth',2)
%                 plot(time*1000,mean(bsxfun(@(x,y) x-y, [tracedata(drugidxes).pre_y], prey0s(drugidxes))')*1000,'r-','LineWidth',2)
                
                plot(time*1000,[tracedata(ctrlidxes).pre_y]'*1000,'k-','Color',[.8 .8 .8])
                plot(time*1000,nanmean([tracedata(ctrlidxes).pre_y]')*1000,'k-','LineWidth',2)



                axis tight
                xlabel('time (ms)')
                ylabel('pre Vm (mV)')
%                 subplot(2,1,2)
                
                set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Position',[1/xcm 1/ycm 1-2/xcm 1-2/ycm],'Xtick',[-valtozok.baselinelength:valtozok.baselinelength:valtozok.psplength]*1000)
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 xsize/dpi ysize/dpi])
                print(gcf,[dirs.figuresdirnow,xlsdata(prenum).ID,'-to-',xlsdata(findpostidxes(potpostnum)).ID,'-ctrl-pre.jpg'],'-djpeg',['-r',num2str(dpi)])
                clf
                hold on
                plot(time*1000,(bsxfun(@(x,y) x-y, [tracedata(ctrlidxes).post_y], y0s(ctrlidxes))')*1000,'k-','Color',[.8 .8 .8])
                plot(time*1000,nanmean(bsxfun(@(x,y) x-y, [tracedata(ctrlidxes).post_y], y0s(ctrlidxes))')*1000,'k-','LineWidth',2)
                axis tight
                ylimits=get(gca,'Ylim');
                maxval=max(max(nanmean(bsxfun(@(x,y) x-y, [tracedata(ctrlidxes).post_y], y0s(ctrlidxes))')*1000));
                minval=min(min(nanmean(bsxfun(@(x,y) x-y, [tracedata(ctrlidxes).post_y], y0s(ctrlidxes))')*1000));
                dval=maxval-minval;
                ylimits(1)=max(ylimits(1),-dval*3);
                ylimits(2)=min(ylimits(2),dval*3);
                set(gca,'Ylim',ylimits);
                xlabel('time (ms)')
                ylabel('post Vm (mV)')
                
                set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Position',[1/xcm 1/ycm 1-2/xcm 1-2/ycm],'Xtick',[-valtozok.baselinelength:valtozok.baselinelength:valtozok.psplength]*1000)
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 xsize/dpi ysize/dpi])
                print(gcf,[dirs.figuresdirnow,xlsdata(prenum).ID,'-to-',xlsdata(findpostidxes(potpostnum)).ID,'-ctrl-post.jpg'],'-djpeg',['-r',num2str(dpi)])
                figure(2)
                clf
%                 subplot(3,1,3)
                hold on
                [~,xout]=hist(y0s([ctrlidxes]));
                [nc,xout]=hist(y0s(ctrlidxes),xout);
                hb=bar(xout*1000,[nc]','grouped');
%                 set(hb(1),'FaceColor',[0 0 0])
%                 set(hb(2),'FaceColor',[1 0 0])
                xlabel('post V0 (mV)')
                ylabel('count')
                set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Position',[1/xcm 1/ycm 1-2/xcm 1-2/ycm])
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 xsize/dpi ysize/dpi])
                print(gcf,[dirs.figuresdirnow,xlsdata(prenum).ID,'-to-',xlsdata(findpostidxes(potpostnum)).ID,'-ctrl-v0hist.jpg'],'-djpeg',['-r',num2str(dpi)])
            % dividing sweeps according to drugs
            for drugnum=1:xlsdata(prenum).drugnum
                %%
                drugname=xlsdata(prenum).drugdata(drugnum).DrugName;
                ctrlidxes=find([tracedata.pre_realtime]<xlsdata(prenum).drugdata(drugnum).DrugWashinTime);
                drugidxes=find([tracedata.pre_realtime]>xlsdata(prenum).drugdata(drugnum).DrugWashinTime+valtozok.drugwashintime);
                y0s=[tracedata.post_y0];
                prey0s=[tracedata.pre_y0];
                y0stds=[tracedata.post_y0sd];
                if length(drugidxes)>3
                   
                    
                % removing noisy sweeps
                zeroed_post_y=bsxfun(@(x,y) x-y, [tracedata.post_y], y0s);
                voltmar=0;
                while voltmar==0 | any(median(ctrldiffs)+3*std(ctrldiffs)<ctrldiffs) | any(median(ctrldiffs)-3*std(ctrldiffs)>ctrldiffs)
                    ctrldiffs=bsxfun(@(x,y) x-y, [zeroed_post_y(:,ctrlidxes)], nanmean(zeroed_post_y(:,ctrlidxes),2));
                    ctrldiffs=max(abs(ctrldiffs));
                    if any(median(ctrldiffs)+3*std(ctrldiffs)<ctrldiffs)
                        [~,idxtodel]=max(ctrldiffs);
                        ctrlidxes(idxtodel)=[];
                    end
                    if any(median(ctrldiffs)-3*std(ctrldiffs)>ctrldiffs)
                        [~,idxtodel]=min(ctrldiffs);
                        ctrlidxes(idxtodel)=[];
                    end
                    voltmar=1;
                end
                voltmar=0;
                while voltmar==0 | any(median(drugdiffs)+3*std(drugdiffs)<drugdiffs) | any(median(drugdiffs)-3*std(drugdiffs)>drugdiffs)
                    drugdiffs=bsxfun(@(x,y) x-y, [zeroed_post_y(:,drugidxes)], nanmean(zeroed_post_y(:,drugidxes),2));
                    drugdiffs=max(abs(drugdiffs));
                    if any(median(drugdiffs)+3*std(drugdiffs)<drugdiffs)
                        [~,idxtodel]=max(drugdiffs);
                        drugidxes(idxtodel)=[];
                    end
                    if any(median(drugdiffs)-3*std(drugdiffs)>drugdiffs)
                        [~,idxtodel]=min(drugdiffs);
                        drugidxes(idxtodel)=[];
                    end
                    voltmar=1;
                end
                % deleting outlying y0 values
                while any(median(y0s(ctrlidxes))+3*std(y0s(ctrlidxes))<y0s(ctrlidxes)) | any(median(y0s(ctrlidxes))-3*std(y0s(ctrlidxes))>y0s(ctrlidxes))
                    if any(median(y0s(ctrlidxes))+3*std(y0s(ctrlidxes))<y0s(ctrlidxes))
                        [~,idxtodel]=max(y0s(ctrlidxes));
                        ctrlidxes(idxtodel)=[];
                    end
                    if any(median(y0s(ctrlidxes))-3*std(y0s(ctrlidxes))>y0s(ctrlidxes))
                        [~,idxtodel]=min(y0s(ctrlidxes));
                        ctrlidxes(idxtodel)=[];
                    end
                end
                
                while any(median(y0s(drugidxes))+3*std(y0s(drugidxes))<y0s(drugidxes)) | any(median(y0s(drugidxes))-3*std(y0s(drugidxes))>y0s(drugidxes))
                    if any(median(y0s(drugidxes))+3*std(y0s(drugidxes))<y0s(drugidxes))
                        [~,idxtodel]=max(y0s(drugidxes));
                        drugidxes(idxtodel)=[];
                    end
                    if any(median(y0s(drugidxes))-3*std(y0s(drugidxes))>y0s(drugidxes))
                        [~,idxtodel]=min(y0s(drugidxes));
                        drugidxes(idxtodel)=[];
                    end
                end
                
                while abs(nanmean(y0s(ctrlidxes))-nanmean(y0s(drugidxes)))>valtozok.maxy0baselinedifference
                    if length(ctrlidxes)> length(drugidxes)
                        if nanmean(y0s(ctrlidxes))>nanmean(y0s(drugidxes))
                            [~,idxtodel]=max(y0s(ctrlidxes));
                            ctrlidxes(idxtodel)=[];
                        else
                            [~,idxtodel]=min(y0s(ctrlidxes));
                            ctrlidxes(idxtodel)=[];
                        end
                    else
                        if nanmean(y0s(drugidxes))>nanmean(y0s(ctrlidxes))
                            [~,idxtodel]=max(y0s(drugidxes));
                            drugidxes(idxtodel)=[];
                        else
                            [~,idxtodel]=min(y0s(drugidxes));
                            drugidxes(idxtodel)=[];
                        end
                    end
%                     disp([mean(y0s(ctrlidxes)),mean(y0s(drugidxes))])
                end
%                 figure(1)
%                 clf
%                 [~,xout]=hist(y0s([ctrlidxes,drugidxes]));
%                 [nc,xout]=hist(y0s(ctrlidxes),xout);
%                 [nd,xout]=hist(y0s(drugidxes),xout);
%                 hb=bar(xout*1000,[nc;nd]','grouped');
%                 set(hb(1),'FaceColor',[0 0 0])
%                 set(hb(2),'FaceColor',[1 0 0])
%                 xlabel('post V0 (mV)')
%                 ylabel('count')
%                 [x,~]=ginput(2);
%                 
%                 
%                  % ctrl and drug v0s have to be in the same range
%                     minY0=min(x)/1000;
%                     maxY0=max(x)/1000;
%                     ctrlidxes(y0s(ctrlidxes)<minY0 | y0s(ctrlidxes)>maxY0)=[];
%                     drugidxes(y0s(drugidxes)<minY0 | y0s(drugidxes)>maxY0)=[];
%                     
                    % ctrl and drug v0s have to be in the same range
                    
                    
                    % pairing ctrl and drug y0 values
                    
                    
                    % pairing ctrl and drug y0 values
                    
                    
                figure(1)
                clf
                
%                 subplot(2,1,1)
                hold on
%                 plot(time*1000,(bsxfun(@(x,y) x-y, [tracedata(ctrlidxes).pre_y], prey0s(ctrlidxes))')*1000,'k-','Color',[.8 .8 .8])
%                 plot(time*1000,(bsxfun(@(x,y) x-y, [tracedata(drugidxes).pre_y], prey0s(drugidxes))')*1000,'r-','Color',[.9 .4 .4])
%                 plot(time*1000,mean(bsxfun(@(x,y) x-y, [tracedata(ctrlidxes).pre_y], prey0s(ctrlidxes))')*1000,'k-','LineWidth',2)
%                 plot(time*1000,mean(bsxfun(@(x,y) x-y, [tracedata(drugidxes).pre_y], prey0s(drugidxes))')*1000,'r-','LineWidth',2)
                
                plot(time*1000,[tracedata(ctrlidxes).pre_y]'*1000,'k-','Color',[.8 .8 .8])
                plot(time*1000,[tracedata(drugidxes).pre_y]'*1000,'r-','Color',[.9 .4 .4])
                plot(time*1000,nanmean([tracedata(ctrlidxes).pre_y]')*1000,'k-','LineWidth',2)
                plot(time*1000,nanmean([tracedata(drugidxes).pre_y]')*1000,'r-','LineWidth',2)



                axis tight
                xlabel('time (ms)')
                ylabel('pre Vm (mV)')
%                 subplot(2,1,2)
                
                set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Position',[1/xcm 1/ycm 1-2/xcm 1-2/ycm],'Xtick',[-valtozok.baselinelength:valtozok.baselinelength:valtozok.psplength]*1000)
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 xsize/dpi ysize/dpi])
                print(gcf,[dirs.figuresdirnow,xlsdata(prenum).ID,'-to-',xlsdata(findpostidxes(potpostnum)).ID,'-',drugname,'-pre.jpg'],'-djpeg',['-r',num2str(dpi)])
                
                clf
                hold on
                plot(time*1000,(bsxfun(@(x,y) x-y, [tracedata(ctrlidxes).post_y], y0s(ctrlidxes))')*1000,'k-','Color',[.8 .8 .8])
                plot(time*1000,(bsxfun(@(x,y) x-y, [tracedata(drugidxes).post_y], y0s(drugidxes))')*1000,'r-','Color',[.9 .4 .4])
                plot(time*1000,nanmean(bsxfun(@(x,y) x-y, [tracedata(ctrlidxes).post_y], y0s(ctrlidxes))')*1000,'k-','LineWidth',2)
                plot(time*1000,nanmean(bsxfun(@(x,y) x-y, [tracedata(drugidxes).post_y], y0s(drugidxes))')*1000,'r-','LineWidth',2)
                axis tight
                ylimits=get(gca,'Ylim');
                maxval=max(max(nanmean(bsxfun(@(x,y) x-y, [tracedata(ctrlidxes).post_y], y0s(ctrlidxes))')*1000),max(nanmean(bsxfun(@(x,y) x-y, [tracedata(drugidxes).post_y], y0s(drugidxes))')*1000));
                minval=min(min(nanmean(bsxfun(@(x,y) x-y, [tracedata(ctrlidxes).post_y], y0s(ctrlidxes))')*1000),min(nanmean(bsxfun(@(x,y) x-y, [tracedata(drugidxes).post_y], y0s(drugidxes))')*1000));
                dval=maxval-minval;
                ylimits(1)=max(ylimits(1),-dval*3);
                ylimits(2)=min(ylimits(2),dval*3);
                set(gca,'Ylim',ylimits);
                xlabel('time (ms)')
                ylabel('post Vm (mV)')
                
                set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Position',[1/xcm 1/ycm 1-2/xcm 1-2/ycm],'Xtick',[-valtozok.baselinelength:valtozok.baselinelength:valtozok.psplength]*1000)
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 xsize/dpi ysize/dpi])
                print(gcf,[dirs.figuresdirnow,xlsdata(prenum).ID,'-to-',xlsdata(findpostidxes(potpostnum)).ID,'-',drugname,'-post.jpg'],'-djpeg',['-r',num2str(dpi)])
                figure(2)
                clf
%                 subplot(3,1,3)
                hold on
                [~,xout]=hist(y0s([ctrlidxes,drugidxes]));
                [nc,xout]=hist(y0s(ctrlidxes),xout);
                [nd,xout]=hist(y0s(drugidxes),xout);
                hb=bar(xout*1000,[nc;nd]','grouped');
                set(hb(1),'FaceColor',[0 0 0])
                set(hb(2),'FaceColor',[1 0 0])
                xlabel('post V0 (mV)')
                ylabel('count')
                set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Position',[1/xcm 1/ycm 1-2/xcm 1-2/ycm])
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 xsize/dpi ysize/dpi])
                print(gcf,[dirs.figuresdirnow,xlsdata(prenum).ID,'-to-',xlsdata(findpostidxes(potpostnum)).ID,'-',drugname,'-v0hist.jpg'],'-djpeg',['-r',num2str(dpi)])
%                 pause
                end
            end
            
            end
        end
    end
    

end
return

%% 
%plotting IV

betumeret=8;
axesvastagsag=2;


xinch=xcm/2.54;
yinch=ycm/2.54;
xsize=dpi*xinch;
ysize=dpi*yinch;

sweeplevonasokeddig=[2,1,1,2,3,4,0,2,0,0,0,1,1,1,1,2,2,8,4,1,1,1,3,2,2,3,2,0,1,1,2,1,3,3,3,2,3,2,1,1,2,4];
sweeplevonasok=NaN(size(xlsdata));
sweeplevonasok(1:length(sweeplevonasokeddig))=sweeplevonasokeddig;
dothesecond=zeros(size(sweeplevonasok));

for prenum=1:length(xlsdata)
    fname=[xlsdata(prenum).ID,'.mat'];
    HEKAfname=xlsdata(prenum).HEKAfname;
    if any(strfind(HEKAfname,','))
        comma=strfind(HEKAfname,',');
        HEKAfname=HEKAfname(1:comma(1)-1);
    end
    load([locations.tgtardir,'MATLABdata/IV/',xlsdata(prenum).setup,'/',HEKAfname,'.mat']);
    
    gsc=xlsdata(prenum).G_S_C;
    hyps=strfind(gsc,'_');
    commas=strfind(gsc,',');
    g=num2str(gsc(1:hyps(1)-1));
    if isempty(commas)
        s=num2str(gsc(hyps(1)+1:hyps(2)-1));
    elseif dothesecond(prenum)==1 & length(commas)==1
        s=num2str(gsc(commas(1)+1:hyps(2)-1));
    elseif dothesecond(prenum)==1 & length(commas)>1
        s=num2str(gsc(commas(1)+1:commas(2)-1));
    else
        s=num2str(gsc(hyps(1)+1:commas(1)-1));
    end
    c=num2str(gsc(hyps(2)+1:end));
    iv=iv.(['g',g,'_s',s,'_c',c]);
    si=mode(diff(iv.time));
    [b,a]=butter(3,15000/(1/mode(diff(iv.time)))/2,'low');
    [bb,aa]=butter(3,1000/(1/mode(diff(iv.time)))/2,'low');
%     for ii=1:iv.sweepnum
%         if ii<5
%             iv.(['v',num2str(ii)])=filtfilt(bb,aa,iv.(['v',num2str(ii)]));
%         else
%             iv.(['v',num2str(ii)])=filtfilt(b,a,iv.(['v',num2str(ii)]));
%         end
%     end
if isnan(sweeplevonasok((prenum)))
    figure(3)
    clf
    subplot(5,1,1)
    hold on;
    plot(iv.time,iv.v1,'k-','LineWidth',2)
    plot(iv.time,iv.(['v',num2str(iv.sweepnum-4)]),'k-','LineWidth',2);
    axis tight
    title(xlsdata(prenum).ID)
    subplot(5,1,2)
    hold on;
    plot(iv.time,iv.v1,'k-','LineWidth',2)
    plot(iv.time,iv.(['v',num2str(iv.sweepnum-3)]),'k-','LineWidth',2);
    axis tight
    title(['Ca: ',num2str(xlsdata(prenum).Ca),' mM    Mg:',num2str(xlsdata(prenum).Mg),' mM'])
    subplot(5,1,3)
    hold on;
    plot(iv.time,iv.v1,'k-','LineWidth',2)
    plot(iv.time,iv.(['v',num2str(iv.sweepnum-2)]),'k-','LineWidth',2);
    axis tight
    subplot(5,1,4)
    hold on;
    plot(iv.time,iv.v1,'k-','LineWidth',2)
    plot(iv.time,iv.(['v',num2str(iv.sweepnum-1)]),'k-','LineWidth',2);
    axis tight
    subplot(5,1,5)
    hold on;
    plot(iv.time,iv.v1,'k-','LineWidth',2)
    plot(iv.time,iv.(['v',num2str(iv.sweepnum)]),'k-','LineWidth',2);
    axis tight
    pause
    sweeplevonasok(prenum)=0;
end

% 
    figure(33)
    clf
    hold on
    plot(iv.time,filtfilt(bb,aa,iv.v1),'k-','LineWidth',2)
    plot(iv.time,filtfilt(b,a,iv.(['v',num2str(iv.sweepnum-(sweeplevonasok(prenum)))])),'k-','LineWidth',2);
    axis tight
    ylim([-.13 .050])
    xlim([0 1])
%     ylimitek(:,i)=get(gca,'Ylim');
    set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Position',[1/xcm 1/ycm 1-2/xcm 1-2/ycm])
    axis off
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 xsize/dpi ysize/dpi])
    print(gcf,[dirs.figuresdir,'/IVs/IV_',xlsdata(prenum).ID,'.jpg'],'-djpeg',['-r',num2str(dpi)])
end
figure(1)
clf
hold on

plot([.5 .5],[-.050 .0],'k-','LineWidth',5)
plot([.5 .6],[-.050 -.050],'k-','LineWidth',5)
ylim([-.13 .050])
    xlim([0 1])
axis off

set(gca,'Position',[0 0 1 1])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 xsize/dpi ysize/dpi])
print(gcf,[dirs.figuresdir,'/IVs/IVscalebar_100ms_50mV.jpg'],'-djpeg',['-r',num2str(dpi)])

% %% making the big picture
% figfiles=dir(dirs.figuresdir);
% figfiles([figfiles.isdir])=[];
% IVfigfiles=dir([dirs.figuresdir,'IVs/']);
% IVfigfiles([IVfigfiles.isdir])=[];
% for xlsnum=1:length(xlsdata)
%     for filenum=1:length(figfiles)
%         if any(strfind(figfiles(filenum).name,xlsdata(xlsnum).ID))
%             if strfind(figfiles(filenum).name,xlsdata(xlsnum).ID)==1
%                 figfiles(filenum).prexlsidx=xlsnum;
%             else
%                 figfiles(filenum).postxlsidx=xlsnum;
%             end
%             if any(strfind(figfiles(filenum).name,'GJ'))
%                 figfiles(filenum).gj=1;
%             else
%                 figfiles(filenum).gj=0;
%             end
%             if any(strfind(figfiles(filenum).name,'pre'))
%                 figfiles(filenum).pre=1;
%                 figfiles(filenum).post=0;
%             elseif any(strfind(figfiles(filenum).name,'post'))
%                 figfiles(filenum).pre=0;
%                 figfiles(filenum).post=1;
%             else
%                 figfiles(filenum).pre=0;
%                 figfiles(filenum).post=0;
%             end
%             if any(strfind(figfiles(filenum).name,'ctrl'))
%                 figfiles(filenum).ctrl=1;
%             else
%                 figfiles(filenum).ctrl=0;
%             end
%             if any(strfind(figfiles(filenum).name,'v0hist'))
%                 figfiles(filenum).v0hist=1;
%             else
%                 figfiles(filenum).v0hist=0;
%             end
%         end
%     end
%     xlsdata(xlsnum).ivfileidx=find(strcmp(['IV_',xlsdata(xlsnum).ID,'.jpg'],{IVfigfiles.name}));
% end
% for xlsnum=1:length(xlsdata)
%        figs=struct;
%  
%     neededfigidxes=find([figfiles.prexlsidx]==xlsnum|[figfiles.postxlsidx]==xlsnum);
%     neededxlsnums=unique([[figfiles(neededfigidxes).prexlsidx],[figfiles(neededfigidxes).postxlsidx]]);
%     neededxlsnums(neededxlsnums==xlsnum)=[];
%     neededxlsnums=[xlsnum,neededxlsnums];
% 
%     for celli=1:length(neededxlsnums)
%         figs(celli).IV=imread([dirs.figuresdir,'IVs/',IVfigfiles(xlsdata(neededxlsnums(celli)).ivfileidx).name]);
%         for precelli=1:length(neededxlsnums)
%             if celli==precelli;
%                 figidx=find([figfiles.ctrl]==1 &[figfiles.gj]==0& [figfiles.pre]==1 & [figfiles.v0hist]==0 &[figfiles.prexlsidx]==neededxlsnums(celli));
%                 if ~isempty(figidx)
%                     figidx=figidx(1);
%                 end
%             else
%                 figidx=find([figfiles.ctrl]==1 &[figfiles.gj]==0& [figfiles.post]==1 & [figfiles.v0hist]==0& [figfiles.prexlsidx]==neededxlsnums(precelli)& [figfiles.postxlsidx]==neededxlsnums(celli));
%             end
%             if ~isempty(figidx)
%                 figs(celli).pp(precelli).fig=imread([dirs.figuresdir,figfiles(figidx).name]);
%             end
%             if celli==precelli;
%                 figidx=find([figfiles.gj]==1& [figfiles.pre]==1& [figfiles.v0hist]==0 &[figfiles.prexlsidx]==neededxlsnums(celli));
%                 if ~isempty(figidx)
%                     figidx=figidx(1);
%                 end
%             else
%                 figidx=find([figfiles.gj]==1& [figfiles.post]==1 & [figfiles.v0hist]==0&[figfiles.prexlsidx]==neededxlsnums(precelli)& [figfiles.postxlsidx]==neededxlsnums(celli));
%             end
%             if ~isempty(figidx)
%                 figs(celli).gj(precelli).fig=imread([dirs.figuresdir,figfiles(figidx).name]);
%             end
%             
%         end
%         
%     end
% 
%     figure(11111)
%     clf
%     tight_subplot(length(figs),length(figs)*2+1,0,0,0);
%     for celli=1:length(neededxlsnums)
%         subplot(length(figs),length(figs)*2+1,(celli-1)*(length(figs)*2+1)+1)
% %         set(gca,'Position',[.1 (celli-1)/(length(figs)) .1 .1])
%         image(figs(celli).IV)
%         box off
%         axis off
%         for ppi=1:length(figs(celli).pp)
%             subplot(length(figs),length(figs)*2+1,(celli-1)*(length(figs)*2+1)+1+ppi)
%             image(figs(celli).pp(ppi).fig)
%             box off
%             axis off
% %             pause
%         end
%         if isfield(figs,'gj')
%             for gji=1:length(figs(celli).gj)
%                 subplot(length(figs),length(figs)*2+1,(celli-1)*(length(figs)*2+1)+1+ppi+gji)
%                 image(figs(celli).gj(gji).fig)
%                 box off
%                 axis off
%                 %             pause
%             end
%         end
%     end
%     pause
% end
