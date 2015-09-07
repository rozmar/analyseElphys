function aE_findevents(valtozok,dirs)
files=dir(dirs.bridgeddir);
files([files.isdir])=[];
for filenum=length(files):-1:1%1:length(files) % végigmegyünk az összes file-n
    a=dir([dirs.eventdir,files(filenum).name]);
    if isempty(a) || valtozok.overwrite==1
        load([dirs.bridgeddir,files(filenum).name]);
        eventdata=struct;
        threshdys=NaN(size(bridgeddata));
        for sweepnum=1:length(bridgeddata) % threshold kalkulálása
            if strcmp(bridgeddata(sweepnum).channellabel(1:end-2),'Vmon') && strcmp(stimdata(sweepnum).Amplifiermode,'C-Clamp')
                sdvalue=aE_calculatebaselineSD(0:bridgeddata(sweepnum).si:(length(bridgeddata(sweepnum).y)-1)*bridgeddata(sweepnum).si,bridgeddata(sweepnum).y,valtozok);
                threshdys(sweepnum)=sdvalue;
            end
            progressbar(sweepnum/length(bridgeddata),[],[],[]);
        end % threshold kalkulálása eddig
        threshdysold=threshdys;
        for sweepnum=1:length(bridgeddata) % küszöbértékek átlagolása, ezzel robosztusabb lesz a módszer
            idxes=find([bridgeddata.realtime]<bridgeddata(sweepnum).realtime+valtozok.threshholdaveragetime & [bridgeddata.realtime]>bridgeddata(sweepnum).realtime-valtozok.threshholdaveragetime);
            threshdys(sweepnum)=nanmean(threshdysold(idxes));
        end
        baselinesds=threshdys/valtozok.eventminsdval;
        
        %
        % itt kezdődik az események keresése
        for sweepnum=1:length(bridgeddata)
            if strcmp(bridgeddata(sweepnum).channellabel(1:4),'Vmon') && strcmp(stimdata(sweepnum).Amplifiermode,'C-Clamp')
                %serkentő események keresése
                
                
                step=round(valtozok.steptime/bridgeddata(sweepnum).si);
                stepback=ceil(step/2);
                diffmovingstep=round(valtozok.diffmovingt/bridgeddata(sweepnum).si);
                filtmovingstep=round(valtozok.filtermovingtime/bridgeddata(sweepnum).si);
                time=[0:bridgeddata(sweepnum).si:(length(bridgeddata(sweepnum).y)-1)*bridgeddata(sweepnum).si];
                y=bridgeddata(sweepnum).y;
                dyraw=diff(y)/mode(diff(time));
                yfilt=imfilter(y,fspecial('average', [1,filtmovingstep]),'symmetric');
                dy=diff(yfilt)/mode(diff(time));
                dyf=imfilter(dy,fspecial('average', [1,diffmovingstep]));
                for tidx=1:length(dyf)-1
                    if dyf(tidx)*dyf(tidx+1)<0
                        dyf(tidx)=0;
                    end
                end
                valtozok.threshdy=threshdys(sweepnum);
                tempampl=abs(dyf)>valtozok.threshdy;
                tempampl=bwareaopen(tempampl,3);
                [amplmask,amplnum]=bwlabel(tempampl);
                
                ittlehetpsp=true(size(yfilt));
                segmentaverages=[];
                seghossz=diff(stimdata(sweepnum).segmenths);
                for segnum=1:length(seghossz)
                    segmentaverages(segnum)=mean(y(stimdata(sweepnum).segmenths(segnum)+round(seghossz(segnum)/2):stimdata(sweepnum).segmenths(segnum+1)-1));
                end
                for segnum=1:length(seghossz)-1
                    if segmentaverages(segnum)>segmentaverages(segnum+1)
                        ittlehetpsp(stimdata(sweepnum).segmenths(segnum+1):find(yfilt(stimdata(sweepnum).segmenths(segnum+1):end)<segmentaverages(segnum+1),1,'first')+stimdata(sweepnum).segmenths(segnum+1))=false;
                    elseif segmentaverages(segnum)<segmentaverages(segnum+1)
                        ittlehetpsp(stimdata(sweepnum).segmenths(segnum+1):find(yfilt(stimdata(sweepnum).segmenths(segnum+1):end)>segmentaverages(segnum+1),1,'first')+stimdata(sweepnum).segmenths(segnum+1))=false;
                    end
                end
                
                
                
                if valtozok.plotit==1
                    figure(11)%
                    clf
                    subplot(2,1,1)
                    hold on
                    plot(time,yfilt,'k-')
                    subplot(2,1,2)
                    plot(time(1:end-1),dy,'k-');
                    hold on
                    plot(time(1:end-1),dyf,'r-','LineWidth',2);
                end
                amplitudes=[];
                halfwidths=[];
                risetimes=[];
                for eventnum=1:amplnum
                    
                    
                    maxh=find(amplmask==eventnum,1,'last');
                    if dyf(maxh)>0
                        eventszorzo=1;
                        type='ep';
                    else
                        eventszorzo=-1;
                        type='ip';
                    end
                    yfilt=yfilt*eventszorzo;
                    y=y*eventszorzo;
                    
                    while maxh<length(y)-step & any(yfilt(maxh)<yfilt(maxh:maxh+step))
                        maxh=maxh+1;
                    end
                    if maxh+step>length(y)
                        [~,maxht]=max(yfilt(maxh:end));
                    else
                        [~,maxht]=max(yfilt(maxh:maxh+step));
                    end
                    maxh=maxh+maxht-1;
                    onseth=find(amplmask==eventnum,1,'first');
                    while onseth>stepback & any(yfilt(onseth)>yfilt(onseth-stepback:onseth))
                        onseth=onseth-1;
                    end
                    endh=maxh;
                    while endh<length(y)-step & any(yfilt(endh)>yfilt(endh:endh+step))
                        endh=endh+1;
                    end
                    amplitude=yfilt(maxh)-yfilt(onseth);
                    halfh1=find(yfilt(onseth:maxh)<yfilt(onseth)+amplitude/2,1,'last')+onseth-1; %%% A félszélesség számolása jelenleg SI pontosságú
                    halfh2=find(yfilt(maxh:end)<yfilt(onseth)+amplitude/2,1,'first')+maxh-1;%%% A félszélesség számolása jelenleg SI pontosságú
                    if isempty(halfh2)
                        halfh2=length(yfilt);
                    end
                    halfwidth=(halfh2-halfh1)*bridgeddata(sweepnum).si; %%% A félszélesség számolása jelenleg SI pontosságú
                    tempwave=(yfilt(onseth:maxh)-yfilt(onseth));
                    h10=find(tempwave>.1*amplitude,1,'first')+onseth-1;
                    h90=find(tempwave>.9*amplitude,1,'first')+onseth-1;
                    tempwave2=(yfilt(maxh:end)-yfilt(onseth));
                    h10end=find(tempwave2<.1*amplitude,1,'first')+maxh-1;
                    h90end=find(tempwave2<.9*amplitude,1,'first')+maxh-1;
                    risetime=(h90-h10)*bridgeddata(sweepnum).si;
                    maxdvalue=max(dyf(amplmask==eventnum));
                    [apmaxv,apmaxh]=max(y(onseth:endh));
                    apmaxh=apmaxh(1)+onseth-1;
                    temph=apmaxh(1);
                    while temph>1 & (y(temph)>(apmaxv-yfilt(onseth))/2+yfilt(onseth) | dyraw(temph)>valtozok.apthreshval)%| dyraw(temph-1)>valtozok.apthreshval)
                        temph=temph-1;
                    end
                    threshh=temph+1;
                    apamplitude=y(apmaxh)-y(threshh);
                    temph=apmaxh(1);
                    while temph<length(y) & (y(temph)>(apmaxv-yfilt(onseth))/2+yfilt(onseth) | dyraw(temph)>valtozok.apthreshval)%| dyraw(temph+1)>valtozok.apthreshval)
                        temph=temph+1;
                    end
                    apendh=temph-1;
                    apendamplitude=y(apmaxh)-y(apendh);
                    apwidth=(apendh-threshh)*bridgeddata(sweepnum).si;
                    
                    ahalfh1=find(y(threshh:apmaxh)<y(threshh)+amplitude/2,1,'last')+threshh-1; %%% A félszélesség számolása jelenleg SI pontosságú
                    ahalfh2=find(y(apmaxh:end)<y(threshh)+amplitude/2,1,'first')+apmaxh-1;%%% A félszélesség számolása jelenleg SI pontosságú
                    if isempty(ahalfh2)
                        ahalfh2=length(y);
                    end
                    aphalfwidth=(ahalfh2-ahalfh1)*bridgeddata(sweepnum).si; %%% A félszélesség számolása jelenleg SI pontosságú
                    
                    tempwave=(y(threshh:apmaxh)-y(threshh));
                    aph10=find(tempwave>.1*apamplitude,1,'first')+threshh-1;
                    aph90=find(tempwave>.9*apamplitude,1,'first')+threshh-1;
                    aprisetime=(h90-h10)*bridgeddata(sweepnum).si;
                    
                    ezegyap=1;%% 
                    ezegypsp=1;
                    
                    
                    if any(ittlehetpsp(onseth:endh)==0)
                        ezegypsp=0;
                    end
                    if amplitude<valtozok.minampl
                        ezegypsp=0;
                    end
                    if amplitude<valtozok.apampl || apendamplitude<valtozok.apampl
                        ezegyap=0;
                    else
                        ezegypsp=0;
                    end
                    if aphalfwidth<valtozok.minaphw
                        ezegyap=0;
                    end
                    if apwidth>valtozok.maxapwidth
                        ezegyap=0;
                    end
                    if eventszorzo<0
                        ezegyap=0;
                    end
                    if (halfh2-maxh)/(maxh-halfh1)>valtozok.maxdecayriseratio%halfwidth>valtozok.maxfilteredhw
                        ezegypsp=0;
                    end
                    if risetime>valtozok.maxrisetime
                        ezegypsp=0;
                    end
                    
                    %                 figure(1)
                    %                 clf
                    %                   plot(time,yfilt);
                    %                 hold on;
                    %                 plot(time(ittlehetpsp),yfilt(ittlehetpsp),'r-')
                    
                    
%                     if isempty(fieldnames(eventdata)) |  ( eventdata(end).maxh<min(maxh,apmaxh) )
                    if  (ezegypsp==1 || ezegyap==1)  %& amplitude/risetime>valtozok.mindvpdt
                        if isempty(fieldnames(eventdata))
                            NEXT=1;
                        else
                            NEXT=length(eventdata)+1;
                        end
                        
                        eventdata(NEXT).sweepnum=sweepnum;
                        eventdata(NEXT).sweeptime=bridgeddata(sweepnum).realtime;
                        eventdata(NEXT).si=bridgeddata(sweepnum).si;
                        if ezegyap==1
                            eventdata(NEXT).onseth=threshh;
                            eventdata(NEXT).maxh=apmaxh;
                            eventdata(NEXT).endh=apendh;
                            eventdata(NEXT).amplitude=apamplitude*eventszorzo;
                            eventdata(NEXT).halfwidth=aphalfwidth;
                            eventdata(NEXT).h10=aph10;
                            eventdata(NEXT).h90=aph90;
                            eventdata(NEXT).risetime=aprisetime;
                            eventdata(NEXT).type='AP';
                            eventdata(NEXT).maxval=y(apmaxh)*eventszorzo;
                            eventdata(NEXT).maxtime=time(apmaxh)+eventdata(NEXT).sweeptime;
                            eventdata(NEXT).onsettime=time(threshh)+eventdata(NEXT).sweeptime;
                            type='AP';
                            % AP után az AHP alaptt nem keresünk PSP-ket
                            % % %                         innenh=onseth;
                            % % %                         eddigh=endh;
                            % % %                         while eddigh<length(y) &  yfilt(eddigh+1) < yfilt(onseth)
                            % % %                             eddigh=eddigh+1;
                            % % %                         end
                            % % %                         ittlehetpsp(innenh:eddigh)=false;
                        else
                            eventdata(NEXT).onseth=onseth;
                            eventdata(NEXT).maxh=maxh;
                            eventdata(NEXT).endh=endh;
                            eventdata(NEXT).amplitude=amplitude*eventszorzo;
                            eventdata(NEXT).halfwidth=halfwidth;
                            eventdata(NEXT).h10=h10;
                            eventdata(NEXT).h90=h90;
                            eventdata(NEXT).risetime=risetime;
                            eventdata(NEXT).type=type;
                            eventdata(NEXT).maxval=yfilt(maxh)*eventszorzo;
                            eventdata(NEXT).maxtime=time(maxh)+eventdata(NEXT).sweeptime;
                            eventdata(NEXT).onsettime=time(onseth)+eventdata(NEXT).sweeptime;
                        end
                        eventdata(NEXT).maxdvalue=maxdvalue;
                        eventdata(NEXT).baselinesd_of_the_sweep=baselinesds(sweepnum);
                        eventdata(NEXT).baselineval=yfilt(onseth)*eventszorzo;
                        
                        eventdata(NEXT).injectedcurrentatpeak=(stimdata(sweepnum).y(maxh));
                        eventdata(NEXT).injectedrelativecurrentatpeak=(stimdata(sweepnum).y(maxh)-stimdata(sweepnum).y(1));
                        eventdata(NEXT).injectedcurrentatonset=(stimdata(sweepnum).y(onseth));
                        eventdata(NEXT).injectedrelativecurrentatonset=(stimdata(sweepnum).y(onseth)-stimdata(sweepnum).y(1));
                        
                        eventdata(NEXT).stimulated=(eventdata(NEXT).injectedrelativecurrentatpeak+eventdata(NEXT).injectedrelativecurrentatonset)>0;
                        
                        tttttempdiffs=stimdata(sweepnum).segmenths-apmaxh;
                        [~,tidx]=min(abs(tttttempdiffs));
                        eventdata(NEXT).maxtimetosquarepulse=tttttempdiffs(tidx)*bridgeddata(sweepnum).si;
                        
                        tttttempdiffs=stimdata(sweepnum).segmenths-onseth;
                        [~,tidx]=min(abs(tttttempdiffs));
                        eventdata(NEXT).onsettimetosquarepulse=tttttempdiffs(tidx)*bridgeddata(sweepnum).si;
                        
                        if valtozok.plotit==1
                            figure(11)
                            subplot(2,1,1)
                            
                            if strcmp(type,'ep')
                                plot(time(onseth:endh),yfilt(onseth:endh)*eventszorzo,'r-','LineWidth',2)
                            elseif strcmp(type,'ip')
                                plot(time(onseth:endh),yfilt(onseth:endh)*eventszorzo,'b-','LineWidth',2)
                            end
                            if ezegyap==1
                                plot(time(onseth:endh),y(onseth:endh)*eventszorzo,'k-');
                                plot(time([threshh,apmaxh]),y([threshh,apmaxh])*eventszorzo,'ro');
                            end
                        end
                    end
%                     end
                    progressbar(sweepnum/length(bridgeddata),[],eventnum/amplnum,[]);
                    yfilt=yfilt*eventszorzo;
                    y=y*eventszorzo;
                end
                
                if valtozok.plotit==1
                    % % %                 figure(2)
                    % % %                 subplot(3,2,2)
                    % % %                 hist(amplitudes,[0:.01:10])
                    % % %                 subplot(3,2,4)
                    % % %                 hist(halfwidths,[0:.1:20])
                    % % %                 subplot(3,2,6)
                    % % %                 hist(abs(amplitudes./risetimes),[0:.1:20])
                    
                    pause
                end
            end
            progressbar(sweepnum/length(bridgeddata),[],[],[]);
        end
        % itt vegzodik az esemenyek keresese
        i=1;
        maradoidxek=1:length(eventdata);
        while i<length(maradoidxek)-1;
            if eventdata(maradoidxek(i)).maxtime==eventdata(maradoidxek(i+1)).maxtime %eventdata(i).onsettime==eventdata(i+1).onsettime
                if eventdata(maradoidxek(i)).onsettime<eventdata(maradoidxek(i+1)).onsettime
                    maradoidxek(i+1)=[];
                else
                    maradoidxek(i)=[];
                    
                end
            else
                i=i+1;
            end
            progressbar(sweepnum/length(bridgeddata),[],[],i/(length(maradoidxek)-1));
        end
        eventdata=eventdata(maradoidxek);
        save([dirs.eventdir,files(filenum).name],'eventdata','valtozok');
        disp([files(filenum).name,' eventfinding done.. ', num2str(length(eventdata)),' events found'])
    else
        disp([files(filenum).name,' eventfinding already done.. skipped'])
    end
    
end
end