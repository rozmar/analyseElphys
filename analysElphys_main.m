%%
% alapadatok
close all
clear all
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
xlsdata=aE_readxls([dirs.basedir,'cb1elphys.xls']);

%% export Raw data from HEKA
for xlsidx=1:length(xlsdata)
    a=dir([dirs.rawexporteddir,xlsdata(xlsidx).ID,'.mat']);
    if isempty(a) | overwrite==1
        rawdata=HEKAexportbytime_main(xlsdata(xlsidx).HEKAfname,xlsdata(xlsidx).setup,xlsdata(xlsidx).Channel,xlsdata(xlsidx).startT,xlsdata(xlsidx).endT);
        save([dirs.rawexporteddir,xlsdata(xlsidx).ID],'rawdata','xlsdata','xlsidx')
        disp([xlsdata(xlsidx).ID,' done'])
    else
        disp([xlsdata(xlsidx).ID,' already done.. skipped'])
    end
end
xlsdataold=xlsdata;

%% generate PGF data, bridge balancing
overwrite=0;
plotRSvalues=0;
plotbridgedsweeps=0;
RSbaselinelength=.00005; %ennyi időt átlagol össze a feszültség megmérésekor
RSrisetime=.00005; %ennyi időt hagy ki az áram injekcióját követően
poolRStime=30; %ebben az idoablakban atlagolja ossze az RS-t a sweep-ek kozott
files=dir(dirs.rawexporteddir);
files([files.isdir])=[];
for fnum=length(files):-1:1%1:length(files)
    fname=files(fnum).name(1:end-4);
%     load([dirs.rawexporteddir,files(fnum).name],'xlsidx');
    xlsidx=find(strcmp({xlsdata.ID},fname));
    if isempty(xlsidx)
        disp(['xls file és filenevek közti összetűzés'])
        pause
    end
    a=dir([dirs.bridgeddir,xlsdata(xlsidx).ID,'.mat']);
    if isempty(a) | overwrite==1
        load([dirs.rawexporteddir,files(fnum).name]);
        xlsdata=xlsdataold;
        % realtime past midnight error
idxtocorrect=find([rawdata.realtime]<rawdata(1).realtime);
for idxi=1:length(idxtocorrect)
rawdata(idxtocorrect(idxi)).realtime=rawdata(idxtocorrect(idxi)).realtime+24*3600;
end
        %     rawdata=temp.rawdata;
        %     xlsdata=temp.xlsdata;
        %     xlsidx=temp.xlsidx;
        stimdata=aE_generatePGF_calculateRS(rawdata,plotRSvalues,RSbaselinelength,RSrisetime,poolRStime);  % generating PGF data, calculating RS
        %bridge balancing
        bridgeddata=struct;
        for sweepnum=1:length(rawdata)
            time=[0:rawdata(sweepnum).si:rawdata(sweepnum).si*(length(stimdata(sweepnum).y)-1)];
            bridgeddata(sweepnum).y=rawdata(sweepnum).y-stimdata(sweepnum).yforbridge*stimdata(sweepnum).RS;
            bridgeddata(sweepnum).si=rawdata(sweepnum).si;
            bridgeddata(sweepnum).realtime=rawdata(sweepnum).realtime;
            bridgeddata(sweepnum).timertime=rawdata(sweepnum).timertime;
            bridgeddata(sweepnum).channellabel=rawdata(sweepnum).channellabel;
            if plotbridgedsweeps==1 & max(diff(stimdata(sweepnum).y))>100*10^-12
                figure(12211)
                clf
                subplot(3,1,1)
                plot(time,rawdata(sweepnum).y)
                hold on
                plot(time(stimdata(sweepnum).segmenths),rawdata(sweepnum).y(stimdata(sweepnum).segmenths),'ro')
                title('raw data')
                subplot(3,1,2)
                plot(time,stimdata(sweepnum).y)
                ylim([-950*10^-12 950*10^-12])
                title('stimulus')
                subplot(3,1,3)
                plot(time,bridgeddata(sweepnum).y)
                title('bridged data')
                pause
            end
        end
        %bridge balancing
        save([dirs.bridgeddir,xlsdata(xlsidx).ID],'stimdata','bridgeddata','RSbaselinelength','RSrisetime','poolRStime','xlsdata','xlsidx','-v7.3')
        disp([xlsdata(xlsidx).ID,' done (bridge_stim)'])
    else
        disp([xlsdata(xlsidx).ID,' already done .. skipped (bridge_stim)'])
    end
end
% return
% % % %% search for NAN data ezek jottek ki eddig: 1410293rm_4_1_3.mat 1211291rm_4_3_4.mat
% % % startdir=dirs.bridgeddir;
% % % files=dir(startdir);
% % % files([files.isdir])=[];
% % % for i=119: length(files)
% % %     fname=files(i).name;
% % %     load([startdir,fname],'bridgeddata');
% % %     nanperc(i)=sum(isnan([bridgeddata.y]))/length([bridgeddata.y]);
% % %     disp([fname,'  -  ',num2str(nanperc(i))]);
% % %     
% % % end
%% finding events
valtozok.overwrite=0;
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
%% 

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
valtozok.baselinelength=0.05; %s
valtozok.psplength=.3; %s
valtozok.filterorder=3;
valtozok.cutofffreq=1500;
valtozok.drugwashintime=120;
valtozok.maxy0baselinedifference=.0005;
valtozok.discardpostsweepswithap=1;

for prenum=17:length(xlsdata) %going throught potential presynaptic cells
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
        neededapidx=find(prediffs>3 & postdiffs<.065 & postdiffs>.055 & postpostdiffs>.1); %criteria to select presynaptic APs
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
                postsweepnum=find(prehypsweeps(prehypsweepnum).realtime==[posttraces.bridgeddata.realtime]);
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
                postsweepidx=find(pretraces.bridgeddata(neededapsweepnums(preapnum)).realtime==[posttraces.bridgeddata.realtime]);
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
            return
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
                    
                if tracedata(sweepnum).apmaxh>stepback    
                    tracedata(sweepnum).pre_y=filter(b,a,tracedata(sweepnum).pre_y)';
                    tracedata(sweepnum).post_y=filter(b,a,tracedata(sweepnum).post_y)';
                    
                    tracedata(sweepnum).pre_y=tracedata(sweepnum).pre_y(tracedata(sweepnum).apmaxh-stepback:tracedata(sweepnum).apmaxh+stepforward);
                    tracedata(sweepnum).post_y=tracedata(sweepnum).post_y(tracedata(sweepnum).apmaxh-stepback:tracedata(sweepnum).apmaxh+stepforward);
                else
                    
                    filtered=filter(b,a,[tracedata(sweepnum).pre_y(end:-1:1),tracedata(sweepnum).pre_y]);
                    tracedata(sweepnum).pre_y=filtered(end-length(tracedata(sweepnum).pre_y):end)';
                    filtered=filter(b,a,[tracedata(sweepnum).post_y(end:-1:1),tracedata(sweepnum).post_y]);
                    tracedata(sweepnum).post_y=filtered(end-length(tracedata(sweepnum).pre_y):end)';
                   
                    
                    tracedata(sweepnum).pre_y=[nan(stepback-tracedata(sweepnum).apmaxh+1,1);tracedata(sweepnum).pre_y(1:tracedata(sweepnum).apmaxh+stepforward)];
                    tracedata(sweepnum).post_y=[nan(stepback-tracedata(sweepnum).apmaxh+1,1);tracedata(sweepnum).post_y(1:tracedata(sweepnum).apmaxh+stepforward)];
                end
                    
                
                    
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
%% 
%plotting IV

betumeret=8;
axesvastagsag=2;


xinch=xcm/2.54;
yinch=ycm/2.54;
xsize=dpi*xinch;
ysize=dpi*yinch;

sweeplevonasokeddig=[2,1,1,2,3,4,0,2,0,0,0,1,1,1,1,2,2,0,4,1,1,1,3];
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
end

sweeplevonasok(prenum)=0;
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