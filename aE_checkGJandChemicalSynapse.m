function aE_checkGJandChemicalSynapse(valtozok,xlsdata,dirs)
% clear valtozok

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
        
        findpostidxes=find(potpostidx);
        for potpostnum=1:length(findpostidxes) %going throught potential postsynaptic cells
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
            
            tracedataGJ=aE_testforGJ(valtozok,dirs,pretraces,posttraces,preevents,postevents); % finding traces to test GJ coupling
            if length(tracedataGJ)>1 % plotting GJ coupling DATA
                neededtraces=1:length(tracedataGJ);
                time=tracedataGJ(neededtraces(1)).time;
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
                set(gca,'LineWidth',valtozok.plot.axesvastagsag,'FontSize',valtozok.plot.betumeret,'Fontname',valtozok.plot.betutipus,'Position',[1/valtozok.plot.xcm 1/valtozok.plot.ycm 1-2/valtozok.plot.xcm 1-2/valtozok.plot.ycm])%,'Xtick',[-valtozok.baselinelength:valtozok.baselinelength:valtozok.psplength]*1000
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 valtozok.plot.xsize/valtozok.plot.dpi valtozok.plot.ysize/valtozok.plot.dpi])
                print(gcf,[dirs.figuresdirnow,xlsdata(prenum).ID,'-to-',xlsdata(findpostidxes(potpostnum)).ID,'-GJ-pre.jpg'],'-djpeg',['-r',num2str(valtozok.plot.dpi)])
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
                set(gca,'LineWidth',valtozok.plot.axesvastagsag,'FontSize',valtozok.plot.betumeret,'Fontname',valtozok.plot.betutipus,'Position',[1/valtozok.plot.xcm 1/valtozok.plot.ycm 1-2/valtozok.plot.xcm 1-2/valtozok.plot.ycm])%,'Xtick',[-valtozok.baselinelength:valtozok.baselinelength:valtozok.psplength]*1000
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 valtozok.plot.xsize/valtozok.plot.dpi valtozok.plot.ysize/valtozok.plot.dpi])
                print(gcf,[dirs.figuresdirnow,xlsdata(prenum).ID,'-to-',xlsdata(findpostidxes(potpostnum)).ID,'-GJ-post.jpg'],'-djpeg',['-r',num2str(valtozok.plot.dpi)])
            end
            
            [tracedata,time]=aE_testforChemicalSynapse(valtozok,dirs,pretraces,posttraces,preevents,postevents); %extracting spike triggered data
            todel=false(size(tracedata));
            if length(tracedata)>1
                if strcmp(valtozok.prerecordingmode,'C-Clamp')
                    todel=todel|~strcmp({tracedata.pre_Amplifiermode},valtozok.prerecordingmode);
                    for tempi=1:length(tracedata)
                        todel(tempi)=todel(tempi)|any(strfind(tracedata(tempi).pre_channellabel,'Imon'));
                    end
                elseif strcmp(valtozok.prerecordingmode,'V-Clamp')
                    todel=todel|~strcmp({tracedata.pre_Amplifiermode},valtozok.prerecordingmode);
                    for tempi=1:length(tracedata)
                        todel(tempi)=todel(tempi)|any(strfind(tracedata(tempi).pre_channellabel,'Vmon'));
                    end
                end
                
                if strcmp(valtozok.postrecordingmode,'C-Clamp')
                    todel=todel|~strcmp({tracedata.post_Amplifiermode},valtozok.postrecordingmode);
                    for tempi=1:length(tracedata)
                        todel(tempi)=todel(tempi)|any(strfind(tracedata(tempi).post_channellabel,'Imon'));
                    end
                elseif strcmp(valtozok.postrecordingmode,'V-Clamp')
                    todel=todel|~strcmp({tracedata.post_Amplifiermode},valtozok.postrecordingmode);
                    for tempi=1:length(tracedata)
                        todel(tempi)=todel(tempi)|any(strfind(tracedata(tempi).post_channellabel,'Vmon'));
                    end
                end
            end
            tracedata(todel)=[];
            if length(tracedata)>1
                
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
                if length(ctrlidxes)>1 %~isempty(ctrlidxes)
                figure(1)
                clf
                
                hold on
                 plot(time*1000,[tracedata(ctrlidxes).pre_y]'*1000,'k-','Color',[.8 .8 .8])
                plot(time*1000,nanmean([tracedata(ctrlidxes).pre_y]')*1000,'k-','LineWidth',2)
                axis tight
                xlabel('time (ms)')
                ylabel('pre Vm (mV)')
                
                set(gca,'LineWidth',valtozok.plot.axesvastagsag,'FontSize',valtozok.plot.betumeret,'Fontname',valtozok.plot.betutipus,'Position',[1/valtozok.plot.xcm 1/valtozok.plot.ycm 1-2/valtozok.plot.xcm 1-2/valtozok.plot.ycm],'Xtick',[-valtozok.baselinelength:valtozok.baselinelength:valtozok.psplength]*1000)
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 valtozok.plot.xsize/valtozok.plot.dpi valtozok.plot.ysize/valtozok.plot.dpi])
                print(gcf,[dirs.figuresdirnow,xlsdata(prenum).ID,'-to-',xlsdata(findpostidxes(potpostnum)).ID,'-ctrl-pre.jpg'],'-djpeg',['-r',num2str(valtozok.plot.dpi)])


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
                
                set(gca,'LineWidth',valtozok.plot.axesvastagsag,'FontSize',valtozok.plot.betumeret,'Fontname',valtozok.plot.betutipus,'Position',[1/valtozok.plot.xcm 1/valtozok.plot.ycm 1-2/valtozok.plot.xcm 1-2/valtozok.plot.ycm],'Xtick',[-valtozok.baselinelength:valtozok.baselinelength:valtozok.psplength]*1000)
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 valtozok.plot.xsize/valtozok.plot.dpi valtozok.plot.ysize/valtozok.plot.dpi])
                print(gcf,[dirs.figuresdirnow,xlsdata(prenum).ID,'-to-',xlsdata(findpostidxes(potpostnum)).ID,'-ctrl-post.jpg'],'-djpeg',['-r',num2str(valtozok.plot.dpi)])
                figure(2)
                clf
                hold on
                [~,xout]=hist(y0s([ctrlidxes]));
                [nc,xout]=hist(y0s(ctrlidxes),xout);
                hb=bar(xout*1000,[nc]','grouped');
                xlabel('post V0 (mV)')
                ylabel('count')
                set(gca,'LineWidth',valtozok.plot.axesvastagsag,'FontSize',valtozok.plot.betumeret,'Fontname',valtozok.plot.betutipus,'Position',[1/valtozok.plot.xcm 1/valtozok.plot.ycm 1-2/valtozok.plot.xcm 1-2/valtozok.plot.ycm])
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 valtozok.plot.xsize/valtozok.plot.dpi valtozok.plot.ysize/valtozok.plot.dpi])
                print(gcf,[dirs.figuresdirnow,xlsdata(prenum).ID,'-to-',xlsdata(findpostidxes(potpostnum)).ID,'-ctrl-v0hist.jpg'],'-djpeg',['-r',num2str(valtozok.plot.dpi)])
                % dividing sweeps according to drugs
                
                for drugnum=1:xlsdata(prenum).drugnum
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
                        
                        if ~isempty(drugidxes)
                        figure(1)
                        clf
                        hold on
                        plot(time*1000,[tracedata(ctrlidxes).pre_y]'*1000,'k-','Color',[.8 .8 .8])
                        plot(time*1000,[tracedata(drugidxes).pre_y]'*1000,'r-','Color',[.9 .4 .4])
                        plot(time*1000,nanmean([tracedata(ctrlidxes).pre_y]')*1000,'k-','LineWidth',2)
                        plot(time*1000,nanmean([tracedata(drugidxes).pre_y]')*1000,'r-','LineWidth',2)
                        axis tight
                        xlabel('time (ms)')
                        ylabel('pre Vm (mV)')
                        set(gca,'LineWidth',valtozok.plot.axesvastagsag,'FontSize',valtozok.plot.betumeret,'Fontname',valtozok.plot.betutipus,'Position',[1/valtozok.plot.xcm 1/valtozok.plot.ycm 1-2/valtozok.plot.xcm 1-2/valtozok.plot.ycm],'Xtick',[-valtozok.baselinelength:valtozok.baselinelength:valtozok.psplength]*1000)
                        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 valtozok.plot.xsize/valtozok.plot.dpi valtozok.plot.ysize/valtozok.plot.dpi])
                        print(gcf,[dirs.figuresdirnow,xlsdata(prenum).ID,'-to-',xlsdata(findpostidxes(potpostnum)).ID,'-',drugname,'-pre.jpg'],'-djpeg',['-r',num2str(valtozok.plot.dpi)])
                        
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
                        set(gca,'LineWidth',valtozok.plot.axesvastagsag,'FontSize',valtozok.plot.betumeret,'Fontname',valtozok.plot.betutipus,'Position',[1/valtozok.plot.xcm 1/valtozok.plot.ycm 1-2/valtozok.plot.xcm 1-2/valtozok.plot.ycm],'Xtick',[-valtozok.baselinelength:valtozok.baselinelength:valtozok.psplength]*1000)
                        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 valtozok.plot.xsize/valtozok.plot.dpi valtozok.plot.ysize/valtozok.plot.dpi])
                        print(gcf,[dirs.figuresdirnow,xlsdata(prenum).ID,'-to-',xlsdata(findpostidxes(potpostnum)).ID,'-',drugname,'-post.jpg'],'-djpeg',['-r',num2str(valtozok.plot.dpi)])
                        figure(2)
                        clf
                        hold on
                        [~,xout]=hist(y0s([ctrlidxes,drugidxes]));
                        [nc,xout]=hist(y0s(ctrlidxes),xout);
                        [nd,xout]=hist(y0s(drugidxes),xout);
                        hb=bar(xout*1000,[nc;nd]','grouped');
                        set(hb(1),'FaceColor',[0 0 0])
                        set(hb(2),'FaceColor',[1 0 0])
                        xlabel('post V0 (mV)')
                        ylabel('count')
                        set(gca,'LineWidth',valtozok.plot.axesvastagsag,'FontSize',valtozok.plot.betumeret,'Fontname',valtozok.plot.betutipus,'Position',[1/valtozok.plot.xcm 1/valtozok.plot.ycm 1-2/valtozok.plot.xcm 1-2/valtozok.plot.ycm])
                        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 valtozok.plot.xsize/valtozok.plot.dpi valtozok.plot.ysize/valtozok.plot.dpi])
                        print(gcf,[dirs.figuresdirnow,xlsdata(prenum).ID,'-to-',xlsdata(findpostidxes(potpostnum)).ID,'-',drugname,'-v0hist.jpg'],'-djpeg',['-r',num2str(valtozok.plot.dpi)])
                        end
                    end
                end
                end 
            end
        end
    end
    
    
end
