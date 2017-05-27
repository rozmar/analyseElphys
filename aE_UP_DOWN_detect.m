function aE_UP_DOWN_detect(valtozok,dirs,xlsdata)    
overwrite=valtozok.overwrite;
segmentlength=valtozok.segmentlength;
APdelwin=[.002,.01];
for filei=1:length({xlsdata.ID})
        a=dir([dirs.statedir,xlsdata(filei).ID,'.mat']);
        if isempty(a) | overwrite==1
        disp(['exporting state transitions from  ', xlsdata(filei).ID])
        %      [Selection,ok] = listdlg('ListString',{xlsdata.ID},'ListSize',[300 600]); % az XLS file alapján kiválasztjuk, hogy melyik file összes mérésén szeretnénk végigmenni
        temp=load([dirs.bridgeddir,xlsdata(filei).ID]);
        load([dirs.eventdir,xlsdata(filei).ID],'eventdata');
        bridgeddata=temp.bridgeddata;
        %%
        statedata.UP=struct;
        statedata.DOWN=struct;
        for sweep=1:length(bridgeddata)
            progressbar(sweep/length(bridgeddata))
            si=bridgeddata(sweep).si;
            time=[1:length(bridgeddata(sweep).y)]*si-si;
            duration=time(end)-time(1);
            if duration>segmentlength
                Y=bridgeddata(sweep).y;
                %%
                yfiltered=Y;
%                 eventstoclear=find(strcmp({eventdata.type},'AP')&[eventdata.sweepnum]==sweep);
%                 if ~isempty(eventstoclear)
%                     steptodel=round(APdelwin/si);
%                     for eventi=1:length(eventstoclear)
%                         maxh=eventdata(eventstoclear(eventi)).maxh;
%                         yfiltered(max(1,round(maxh-steptodel(1))):min(length(yfiltered),round(maxh+steptodel(2))))=NaN;
%                     end
%                     ittanan=1;
%                     while ~isempty(ittanan)
%                         ittanan=find(isnan(yfiltered),1,'first');
%                         if ~isempty(ittanan)
%                             ezmarnemnan=find(~isnan(yfiltered(ittanan:end)),1,'first')+ittanan;
%                             if ittanan>1
%                                 kezdo=yfiltered(ittanan-1);
%                             else
%                                 kezdo=yfiltered(ezmarnemnan);
%                             end
%                             if ~isempty(ezmarnemnan)
%                                 vegzo=yfiltered(ezmarnemnan);
%                             else
%                                 ezmarnemnan=lengt(yfiltered)+1;
%                                 vegzo=kezdo;
%                             end
%                             idxdifi=ezmarnemnan-ittanan;
%                             vektoride=([1:idxdifi]/idxdifi)*(vegzo-kezdo)+kezdo;
%                             yfiltered(ittanan:ezmarnemnan-1)=vektoride;
%                         end
%                     end
%                 end
                %%
                downsampletimes=round(.0001/si);
                [b,a]=butter(1,1000/(1/si)/2,'low');
                yfiltered=filtfilt(b,a,yfiltered);
                yfiltered=downsample(yfiltered,downsampletimes);
                timefiltered=downsample(time,downsampletimes);
                yfiltered=medfilt1(yfiltered,10);
                movingwindowstep=.1;
                transitionvoltage=NaN(round(1/movingwindowstep),length(yfiltered));
                transitionvoltage3SD=NaN(size(transitionvoltage));
                for segmenti=1:movingwindowstep:round(duration/segmentlength)
                    sor=int8(mod(segmenti*10,1/movingwindowstep)+1);
                    idx=timefiltered>(segmenti-1)*segmentlength & timefiltered<segmenti*segmentlength;
                    y=yfiltered(idx);
                    yhist=sort(y);
                    yhist=yhist(1:round(.99*length(yhist)));
                    medy=median(y);
                    [n,xout]=hist(yhist,[-20:.2:20]/1000+medy);
                    nmax=moving(n,5,'max');
                    nmin=moving(n,5,'min');
                    maxidx=find(diff(nmax)<0,1,'first');
                    maxidx=find(n==nmax(maxidx-1),1,'first');
                    minidx=find(diff(nmin(maxidx+5:end))>0,1,'first')+maxidx+5;
                    if isempty(minidx)
                        minidx=find(diff(nmin(maxidx+5:end))==0,1,'first')+maxidx+5;
                    end
                    minidx=find(n==nmin(minidx-1),1,'first');
                    maxidx2=find(diff(nmax(minidx+5:end))<0,1,'first')+minidx+5;
                    if isempty(maxidx2)
                        maxidx2=find(diff(nmax(minidx+5:end))==0,1,'first')+minidx+5;
                    end
                    maxidx2=find(n(minidx:end)==nmax(maxidx2-1),1,'first')+minidx;
                    [~,minidx]=min(n(maxidx:maxidx2));
                    minidx=minidx+maxidx-1;
                    [~,maxidxnew]=max(n(maxidx:minidx));
                    maxidx=maxidx+maxidxnew-1;
                    %                 new=ones(size( transitionvoltage(idx)))*xout(minidx);
                    %                 transitionvoltage(idx)= nanmean([transitionvoltage(idx),new]);%xout(minidx);
                    % baseline SD method
                    belowbase=yhist(yhist<xout(maxidx));
                    overbase=-(belowbase-2*max(belowbase));
                    base=[belowbase,overbase];
                    baseSD=std(base);
                    TransitionVoltage=xout(maxidx)+baseSD*3;
                    new=ones(1,sum(idx))*TransitionVoltage;
                    newSD=ones( 1,sum(idx))*baseSD*2;
                    transitionvoltage(sor,idx)= new;%xout(minidx);
                    transitionvoltage3SD(sor,idx)= newSD;%xout(minidx);
                    %
                    
%                                     figure(31)
%                                     clf
%                                     subplot(2,1,1)
%                                     plot(y)
%                                     subplot(2,1,2)
%                                     hold on
%                                     plot( xout,n)
%                                     plot( xout(maxidx),n(maxidx),'ro')
%                                     plot( xout(maxidx2),n(maxidx2),'ro')
%                                     plot( xout(minidx),n(minidx),'bo')
%                                     plot( TransitionVoltage,0,'bx')
%                                     hold on
%                                     plot( xout,nmax,'r-')
%                                     pause
                end
                transitionvoltage=nanmean(transitionvoltage);
                transitionvoltage3SD=nanmean(transitionvoltage3SD);
                %% state transition criteria
                %voltage should exceed transition voltage with 1 mV for state
                %change
                state=zeros(size(transitionvoltage));
                putativestatechanges=find(diff(yfiltered>transitionvoltage)==1 | diff(yfiltered<transitionvoltage)==1);
                putativeUPtransitions=find(diff(yfiltered>transitionvoltage+.001)==1);
                putativeDOWNtransitions=find(diff(yfiltered<transitionvoltage-.001)==1);
                startidx=0;
                endidx=0;
                stepback=10/si;
                for upi=1:length(putativeUPtransitions)
                    if putativeUPtransitions(upi)>endidx & putativeUPtransitions(upi)>500
                        startidxx=find(putativestatechanges<putativeUPtransitions(upi),1,'last');
                        startidx=putativestatechanges(startidxx);
                        
                        nowidx=startidx;
                        downidx=find(putativeDOWNtransitions>nowidx,1,'first');
                        while ~isempty(downidx)& (any(find(putativestatechanges>putativeDOWNtransitions(downidx)&putativestatechanges<putativeDOWNtransitions(downidx)+.1/si/downsampletimes))) % the state transition to the other end should last for at least 100 ms
                                nowidx=putativeDOWNtransitions(downidx);
                                downidx=find(putativeDOWNtransitions>nowidx,1,'first');
                        end
                        
                        if ~isempty(downidx)
                            endidxx=find(putativestatechanges<putativeDOWNtransitions(downidx),1,'last');
                            endidx=putativestatechanges(endidxx);
                        else
                            endidx=length(transitionvoltage);
                        end
                        if (endidx-startidx)*si*downsampletimes>200/1000
                            state(startidx:endidx)=1;
                        end
                    end
                end
                startidx=0;
                endidx=0;
                for downi=1:length(putativeDOWNtransitions)
                    if putativeDOWNtransitions(downi)>endidx & putativeDOWNtransitions(downi)>500
                        startidxx=find(putativestatechanges<putativeDOWNtransitions(downi),1,'last');
                        startidx=putativestatechanges(startidxx);

                        nowidx=startidx;
                        upidx=find(putativeUPtransitions>nowidx,1,'first');
                        while ~isempty(upidx)& (any(find(putativestatechanges>putativeUPtransitions(upidx)&putativestatechanges<putativeUPtransitions(upidx)+.1/si/downsampletimes))) % the state transition to the other end should last for at least 100 ms
                                nowidx=putativeUPtransitions(upidx);
                                upidx=find(putativeUPtransitions>nowidx,1,'first');
                        end
                        if ~isempty(upidx)
                            endidxx=find(putativestatechanges<putativeUPtransitions(upidx),1,'last');
                            endidx=putativestatechanges(endidxx);
                        else
                            endidx=length(transitionvoltage);
                        end
                        if (endidx-startidx)*si*downsampletimes>100/1000
                            state(startidx:endidx)=-1;
                        end
                    end
                end

                UPidx=find(state==1);
                [UPl,num]=bwlabel(state==1);
                for upi=1:num
                    if isempty(fieldnames(statedata.UP))
                        NEXT=1;
                    else
                        NEXT=length(statedata.UP)+1;
                    end
                    idx=find(UPl==upi);
                    statedata.UP(NEXT).type='UP';
                    statedata.UP(NEXT).onseth=idx(1)*downsampletimes;
                    statedata.UP(NEXT).endh=(idx(end)-1)*downsampletimes;
                    statedata.UP(NEXT).onsett=time(idx(1)*downsampletimes)+bridgeddata(sweep).realtime;
                    statedata.UP(NEXT).endt=time((idx(end)-1)*downsampletimes)+bridgeddata(sweep).realtime;
                    statedata.UP(NEXT).midt=median(timefiltered(idx));
                    statedata.UP(NEXT).duration=length(idx)*si*downsampletimes;
                    statedata.UP(NEXT).value=mean(yfiltered(idx));
                    statedata.UP(NEXT).sweepnum=sweep;
                end
                
                DOWNidx=find(state==-1);
                [DOWNl,num]=bwlabel(state==-1);
                for downi=1:num
                    if isempty(fieldnames(statedata.DOWN))
                        NEXT=1;
                    else
                        NEXT=length(statedata.DOWN)+1;
                    end
                    idx=find(DOWNl==downi);
                    statedata.DOWN(NEXT).type='DOWN';
                    statedata.DOWN(NEXT).onseth=idx(1)*downsampletimes;
                    statedata.DOWN(NEXT).endh=(idx(end)-1)*downsampletimes;
                    statedata.DOWN(NEXT).onsett=time(idx(1)*downsampletimes)+bridgeddata(sweep).realtime;
                    statedata.DOWN(NEXT).endt=time((idx(end)-1)*downsampletimes)+bridgeddata(sweep).realtime;
                    statedata.DOWN(NEXT).midt=median(timefiltered(idx));
                    statedata.DOWN(NEXT).duration=length(idx)*si*downsampletimes;
                    statedata.DOWN(NEXT).value=mean(yfiltered(idx));
                    statedata.DOWN(NEXT).sweepnum=sweep;
                end
                if ~isempty(fieldnames(statedata.UP)) & ~isempty(fieldnames(statedata.DOWN))
                    for upi=1:length(statedata.UP)
                        [~,idx]=min(abs([statedata.DOWN.midt]-statedata.UP(upi).midt));
                        statedata.UP(upi).dv=statedata.UP(upi).value-statedata.DOWN(idx).value;
                    end
                    for downi=1:length(statedata.DOWN)
                        [~,idx]=min(abs([statedata.UP.midt]-statedata.DOWN(downi).midt));
                        statedata.DOWN(downi).dv=statedata.DOWN(downi).value-statedata.UP(idx).value;
                    end
                end
%                 figure(3)
%                 clf
%                 plot(time,Y,'k-')
%                 hold on
% %                 plot(timefiltered,yfiltered,'b-')
%                 plot(timefiltered(UPidx),yfiltered(UPidx),'r-')
%                 plot(timefiltered(DOWNidx),yfiltered(DOWNidx),'b-')
%                 pause
            end
        end
        save([dirs.statedir,xlsdata(filei).ID],'statedata');
        end
    end
   