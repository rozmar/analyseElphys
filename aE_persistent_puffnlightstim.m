%%  cutting out spikes
for i=1:length(xlsdata)
    xlsdata(i).newID=[xlsdata(i).ID,'_',xlsdata(i).Stimulation];
end
[Selection,ok] = listdlg('ListString',{xlsdata.newID},'ListSize',[300 600]); % az XLS file alapján kiválasztjuk, hogy melyik file összes mérésén szeretnénk végigmenni
% Selection=1;
APdata=struct;
timebefore=.002;
timeafter=.004;
movingaverage=2;

stimonsetvar.baseline=.005;
stimonsetvar.effectlength=.02;
stimonsetvar.dofiltering=1;
stimonsetvar.filtercutoff=10000;
stimonsetvar.filterorder=3;
stimonsetvar.onsetminSD=10;

for xlsnum=1:length(Selection) %going throught potential presynaptic cells
    prenum=Selection(xlsnum);
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
        if isfield(pretraces.bridgeddata,'stimulation') & pretraces.bridgeddata(sweepnum).stimulation(maxh)==1
            extrastimulated=1;
        else
            
            extrastimulated=0;
        end
        
        si=pretraces.bridgeddata(sweepnum).si;
        stepback=round(timebefore/si);
        stepforward=round(timeafter/si);
        if maxh-stepback*2>0 & maxh+stepforward*2<length(pretraces.bridgeddata(sweepnum).y)
            apwave=moving(pretraces.bridgeddata(sweepnum).y(maxh-stepback:maxh),movingaverage);
            dapwave=diff(apwave)/si;
            [~,dvmaxh]=max(dapwave);
            threshh=find(dapwave(1:dvmaxh)<10,1,'last');
            if isempty(threshh)
                threshh=1;
            end
            dddapwave=(diff(apwave(threshh:end),2));
            [~,maxdddhely]=max(dddapwave);
            
            difi=abs(threshh-stepback);
            difi=0;
            maxh=maxh-difi;
            if ~isempty(maxh)
                APdata(APnumi).y=pretraces.bridgeddata(sweepnum).y(maxh-stepback:maxh+stepforward)';
                APdata(APnumi).x=[-stepback*si:si:stepforward*si]';
                APdata(APnumi).baselineval=mean(pretraces.bridgeddata(sweepnum).y(1:round(.01/si)));
                APdata(APnumi).sweepnum=sweepnum;
                APdata(APnumi).starth=maxh-stepback;
                APdata(APnumi).endh=maxh+stepforward;
                APdata(APnumi).maxh=maxh;
                APdata(APnumi).maxtime=maxh*si+pretraces.bridgeddata(sweepnum).realtime;
                APdata(APnumi).threshval=APdata(APnumi).y(threshh);
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
                
                APdata(APnumi).extrastimulated=extrastimulated;
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
APdataold=APdata;

% cutting out stimulation onsets
LIGHTstimdata=struct;
if any(strfind(xlsdata(prenum).Stimulation,'light'))
for sweepnum=1:length(pretraces.bridgeddata)
    if strcmp(pretraces.bridgeddata(sweepnum).channellabel(1:4),'Vmon') & strcmp(pretraces.stimdata(sweepnum).Amplifiermode,'C-Clamp')
        if max(pretraces.bridgeddata(sweepnum).stimulation)>0
            si=pretraces.bridgeddata(sweepnum).si;
            y=pretraces.bridgeddata(sweepnum).y;
            if stimonsetvar.dofiltering==1
                [b,a]=butter(stimonsetvar.filterorder,stimonsetvar.filtercutoff/(1/si)/2,'low');
                y=filtfilt(b,a,y);
            end
            stepback=round(stimonsetvar.baseline/si);
            stepforward=round(stimonsetvar.effectlength/si);
            time=[-stepback*si:si:stepforward*si];
            [stimbw,stimnum]=bwlabel(pretraces.bridgeddata(sweepnum).stimulation);
            for stimi=1:stimnum
                stimstarth=find(stimbw==stimi,1,'first');
                baselineval=mean(y(stimstarth-stepback:stimstarth));
                baselinesd=std(y(stimstarth-stepback:stimstarth));
                onseth1=find(y(stimstarth:end)>baselineval+baselinesd*stimonsetvar.onsetminSD,1,'first');
                onseth2=find(y(stimstarth:end)<baselineval-baselinesd*stimonsetvar.onsetminSD,1,'first');
                if isempty(onseth1)
                    onseth1=nan;
                end
                if isempty(onseth2)
                    onseth2=nan;
                end
                
                onseth=nanmin(onseth1,onseth2);
%                 figure(1)
%                 clf
%                 plot(time,y(stimstarth-stepback:stimstarth+stepforward));
%                 hold on
%                 plot(time(onseth+stepback+1),y(onseth+stimstarth),'ro')
%                 pause
                if isempty(fieldnames(LIGHTstimdata))
                    LIGHTstimdatai=1;
                else
                    LIGHTstimdatai=length(LIGHTstimdata)+1;
                end
                LIGHTstimdata(LIGHTstimdatai).x=time';
                LIGHTstimdata(LIGHTstimdatai).y=y(stimstarth-stepback:stimstarth+stepforward)';
                LIGHTstimdata(LIGHTstimdatai).si=si;
                LIGHTstimdata(LIGHTstimdatai).stimonseth=stimstarth;
                LIGHTstimdata(LIGHTstimdatai).realonseth=stimstarth+onseth;
                LIGHTstimdata(LIGHTstimdatai).onseth=onseth;
                LIGHTstimdata(LIGHTstimdatai).onsettime=onseth*si;
                if isnan(onseth)
                    LIGHTstimdata(LIGHTstimdatai).onsetval=NaN;
                else
                    LIGHTstimdata(LIGHTstimdatai).onsetval=y(stimstarth+onseth);
                end
                LIGHTstimdata(LIGHTstimdatai).sweepnum=sweepnum;
                LIGHTstimdata(LIGHTstimdatai).baselineval=baselineval;
                LIGHTstimdata(LIGHTstimdatai).baselinesd=baselinesd;
                LIGHTstimdata(LIGHTstimdatai).stimonsettime=pretraces.bridgeddata(sweepnum).realtime+stimstarth*si;
            end
        end
    end
end
end
LIGHTstimdataold=LIGHTstimdata;
%% inspecting spikes and onsets
groupbythresh=1;
% group1=[APdataold.extrastimulated]==0;
% group2=[APdataold.extrastimulated]==1;
figure(4)
clf
hold on
plot([APdataold.baselineval],[APdataold.threshval],'ko')
xlabel('holding potential (V)')
ylabel('threshold potential (V)')

hfig3=figure(3);
clf
hold on
for sweepnum=1:length(pretraces.bridgeddata)
    if strcmp(pretraces.bridgeddata(sweepnum).channellabel(1:4),'Vmon')
        idxesnow=find([APdata.sweepnum]==sweepnum);
        %     if ~isempty(idxesnow)
        si=pretraces.bridgeddata(sweepnum).si;
        hossz=length(pretraces.bridgeddata(sweepnum).y);
        time=[0:si:((hossz-1)*si)]+pretraces.bridgeddata(sweepnum).realtime;
        plot(time,pretraces.bridgeddata(sweepnum).y,'k-','LineWidth',2)
        
        %     end
        if isfield(pretraces.bridgeddata,'stimulation') 
            plot(time(pretraces.bridgeddata(sweepnum).stimulation==1),pretraces.bridgeddata(sweepnum).stimulation(pretraces.bridgeddata(sweepnum).stimulation==1)*min(pretraces.bridgeddata(sweepnum).y),'r-','LineWidth',6)
        end
    end
end
if any(strfind(xlsdata(prenum).Stimulation,'puff')) & isfield(pretraces.bridgeddata,'stimulation')
    hfig5=figure(5);
    clf
    hold on
    for sweepnum=1:length(pretraces.bridgeddata)
        if strcmp(pretraces.bridgeddata(sweepnum).channellabel(1:4),'Vmon')
            idxesnow=find([APdata.sweepnum]==sweepnum);
            %     if ~isempty(idxesnow)
            si=pretraces.bridgeddata(sweepnum).si;
            hossz=length(pretraces.bridgeddata(sweepnum).y);
            time=[0:si:((hossz-1)*si)]+pretraces.bridgeddata(sweepnum).realtime;
                plot(time,pretraces.bridgeddata(sweepnum).stimulation,'k-','LineWidth',2)
        end
    end
    hfig5=gca;
    figure(3)
    hfig3=gca;
    linkaxes([hfig3,hfig5],'x')
end
y0=max([pretraces.bridgeddata.y]);
figure(3)
for i=1:length(xlsdata(prenum).drugdata)
    y0=y0+.005;
    startt=xlsdata(prenum).drugdata(i).DrugWashinTime;
    endt=xlsdata(prenum).drugdata(i).DrugWashoutTime;
    if startt<(pretraces.bridgeddata(1).realtime);
        startt=pretraces.bridgeddata(1).realtime;
    end
    if isempty(endt)
        endt=time(end);
    end
    
    plot([startt, endt],[y0,y0],'k-','LineWidth', 10)
    text(startt,y0,xlsdata(prenum).drugdata(i).DrugName,'Color',[1 1 1])
end
title(xlsdata(prenum).ID)
pause
[x,~]=ginput(4);
x=sort(x);
group1=[APdataold.maxtime]>=x(1) & [APdataold.maxtime]<=x(2);
group2=[APdataold.maxtime]>=x(3) & [APdataold.maxtime]<=x(4);
if any(strfind(xlsdata(prenum).Stimulation,'light'))
    stimgroup1=[LIGHTstimdataold.stimonsettime]>=x(1) & [LIGHTstimdataold.stimonsettime]<=x(2);
    stimgroup2=[LIGHTstimdataold.stimonsettime]>=x(3) & [LIGHTstimdataold.stimonsettime]<=x(4);
end
if groupbythresh==1
    figure(4)
    [~,y]=ginput(1);
    group1=[APdataold.threshval]>=y;
    group2=[APdataold.threshval]<=y;
end

APdata=APdataold(group1);
APdata=APdata([APdata.si]==mode([APdata.si]));
figure(2)
clf
h1=subplot(2,4,1);
plot([APdata.x],[APdata.y])
hold on
plot(mean([APdata.x],2),mean([APdata.y],2),'b-','LineWidth',3)
axis tight
h2=subplot(2,4,2);
plot([APdata.dVperdtT],[APdata.dVperdt])
hold on
plot(mean([APdata.dVperdtT],2),mean([APdata.dVperdt],2),'b-','LineWidth',3)
axis tight
h3=subplot(2,4,3);
plot([APdata.ddVperdtT],[APdata.ddVperdt])
hold on
plot(mean([APdata.ddVperdtT],2),mean([APdata.ddVperdt],2),'b-','LineWidth',3)
axis tight
h4=subplot(2,4,4);
plot([APdata.dVperdtV], [APdata.dVperdt])
hold on
plot(mean([APdata.dVperdtV],2),mean([APdata.dVperdt],2),'b-','LineWidth',3)
axis tight
if any(strfind(xlsdata(prenum).Stimulation,'light'))
    LIGHTstimdata=LIGHTstimdataold(stimgroup1);
    LIGHTstimdata=LIGHTstimdata([LIGHTstimdata.si]==mode([LIGHTstimdata.si]));
    figure(1)
    clf
    hh1=subplot(2,2,1);
    plot([LIGHTstimdata.x],bsxfun(@(x,y) x-y,[LIGHTstimdata.y],[LIGHTstimdata.baselineval]))
    hold on
    plot([LIGHTstimdata.onsettime],[LIGHTstimdata.onsetval]-[LIGHTstimdata.baselineval],'bo')
    plot(mean([LIGHTstimdata.x],2),mean(bsxfun(@(x,y) x-y,[LIGHTstimdata.y],[LIGHTstimdata.baselineval]),2),'b-','LineWidth',3)
    hh2=subplot(2,2,2);
    hist([LIGHTstimdata.onsettime],[0:.0001:.005]);
end
figure(3)
% clf
hold on
for sweepnum=1:length(pretraces.bridgeddata)
    if strcmp(pretraces.bridgeddata(sweepnum).channellabel(1:4),'Vmon')
        idxesnow=find([APdata.sweepnum]==sweepnum);
        %     if ~isempty(idxesnow)
        si=pretraces.bridgeddata(sweepnum).si;
        hossz=length(pretraces.bridgeddata(sweepnum).y);
        time=[0:si:((hossz-1)*si)]+pretraces.bridgeddata(sweepnum).realtime;
        %         plot(time,pretraces.bridgeddata(sweepnum).y,'k-')
        for i=1:length(idxesnow)
            hs=[APdata(idxesnow(i)).starth:APdata(idxesnow(i)).endh];
            plot(time(hs),pretraces.bridgeddata(sweepnum).y(hs),'b-','LineWidth',2)
        end
        
        %     end
    end
end

figure(2)

APdata=APdataold(group2);
APdata=APdata([APdata.si]==mode([APdata.si]));
h5=subplot(2,4,5);
plot([APdata.x],[APdata.y])
hold on
plot(mean([APdata.x],2),mean([APdata.y],2),'r-','LineWidth',3)
axis tight
h6=subplot(2,4,6);
plot([APdata.dVperdtT],[APdata.dVperdt])
hold on
plot(mean([APdata.dVperdtT],2),mean([APdata.dVperdt],2),'r-','LineWidth',3)
axis tight
h7=subplot(2,4,7);
plot([APdata.ddVperdtT],[APdata.ddVperdt])
hold on
plot(mean([APdata.ddVperdtT],2),mean([APdata.ddVperdt],2),'r-','LineWidth',3)
axis tight
h8=subplot(2,4,8);
plot([APdata.dVperdtV], [APdata.dVperdt])
hold on
plot(mean([APdata.dVperdtV],2),mean([APdata.dVperdt],2),'r-','LineWidth',3)
axis tight
linkaxes([h1,h5],'xy')
linkaxes([h2,h6],'xy')
linkaxes([h3,h7],'xy')
linkaxes([h4,h8],'xy')
if any(strfind(xlsdata(prenum).Stimulation,'light'))
    LIGHTstimdata=LIGHTstimdataold(stimgroup2);
    LIGHTstimdata=LIGHTstimdata([LIGHTstimdata.si]==mode([LIGHTstimdata.si]));
    figure(1)
    hh3=subplot(2,2,3);
    plot([LIGHTstimdata.x],bsxfun(@(x,y) x-y,[LIGHTstimdata.y],[LIGHTstimdata.baselineval]))
    hold on
    plot([LIGHTstimdata.onsettime],[LIGHTstimdata.onsetval]-[LIGHTstimdata.baselineval],'ro')
    plot(mean([LIGHTstimdata.x],2),mean(bsxfun(@(x,y) x-y,[LIGHTstimdata.y],[LIGHTstimdata.baselineval]),2),'r-','LineWidth',3)
    hh4=subplot(2,2,4);
    hist([LIGHTstimdata.onsettime],[0:.0001:.005]);
    linkaxes([hh1,hh3],'xy')
    linkaxes([hh2,hh4],'xy')
end
%
figure(3)
% clf
hold on
for sweepnum=1:length(pretraces.bridgeddata)
    if strcmp(pretraces.bridgeddata(sweepnum).channellabel(1:4),'Vmon')
        idxesnow=find([APdata.sweepnum]==sweepnum);
        %     if ~isempty(idxesnow)
        si=pretraces.bridgeddata(sweepnum).si;
        hossz=length(pretraces.bridgeddata(sweepnum).y);
        time=[0:si:((hossz-1)*si)]+pretraces.bridgeddata(sweepnum).realtime;
        %         plot(time,pretraces.bridgeddata(sweepnum).y,'k-')
        for i=1:length(idxesnow)
            hs=[APdata(idxesnow(i)).starth:APdata(idxesnow(i)).endh];
            plot(time(hs),pretraces.bridgeddata(sweepnum).y(hs),'r-','LineWidth',2)
        end
        
        plot(time(pretraces.bridgeddata(sweepnum).stimulation==1),pretraces.bridgeddata(sweepnum).stimulation(pretraces.bridgeddata(sweepnum).stimulation==1)*min(pretraces.bridgeddata(sweepnum).y),'r-','LineWidth',6)
        %     end
    end
end

figure(4)
clf
hold on
plot([APdataold.baselineval],[APdataold.threshval],'ko')
plot([APdataold(group1).baselineval],[APdataold(group1).threshval],'bo','MarkerFaceColor',[0 0 1])
plot([APdataold(group2).baselineval],[APdataold(group2).threshval],'ro','MarkerFaceColor',[1 0 0])
xlabel('holding potential (V)')
ylabel('threshold potential (V)')

if any(strfind(xlsdata(prenum).Stimulation,'puff'))
    linkaxes([hfig3,hfig5],'x')
end
