function [stimdata]=aE_generatePGF_calculateRS(rawdata,plotRSvalues,RSbaselinelength,RSrisetime,poolRStime)


%creating PGF data, calculating RS
stimdata=struct;
RSs=[];
RStimes=[];
Ralls=[];
for sweepnum=1:length(rawdata)
    %some sweeps are exported in nAmps and some in Amps..
    if max(abs(rawdata(sweepnum).segmentamplitudes))>10^-8
        rawdata(sweepnum).segmentamplitudes=rawdata(sweepnum).segmentamplitudes/1000000000;
    end
    RS=NaN;
    RStime=NaN;
    Rall=NaN;
    RSbaselinestep=ceil(RSbaselinelength/rawdata(sweepnum).si);
    RSrisestep=ceil(RSrisetime/rawdata(sweepnum).si);
    
    stimdata(sweepnum).y=zeros(size(rawdata(sweepnum).y));
    preamplnum=str2num(rawdata(sweepnum).channellabel(end));
    stimdata(sweepnum).preamplnum=preamplnum;
    stimdata(sweepnum).Amplifiermode=char(rawdata(sweepnum).AmplifierMode(preamplnum));
    time=[0:rawdata(sweepnum).si:rawdata(sweepnum).si*(length(stimdata(sweepnum).y)-1)];
    if any(isnan(rawdata(sweepnum).segmenttimes))
        ennyisegvan=find(isnan(rawdata(sweepnum).segmenttimes),1,'first')-1;
    else
        ennyisegvan=length(rawdata(sweepnum).segmenttimes);
    end
    for segnum=1:ennyisegvan
        stimdata(sweepnum).segmenths(segnum)=find(time<=rawdata(sweepnum).segmenttimes(segnum),1,'last');
    end
    stimdata(sweepnum).segmenths(segnum+1)=length(time);
    for segnum=1:ennyisegvan
        stimdata(sweepnum).y(stimdata(sweepnum).segmenths(segnum):stimdata(sweepnum).segmenths(segnum+1))=rawdata(sweepnum).segmentamplitudes(segnum);
    end
    if min(stimdata(sweepnum).y)==0 && max(stimdata(sweepnum).y)==0
        stimdata(sweepnum).y=stimdata(sweepnum).y+rawdata(sweepnum).Amplifierholding(preamplnum);
    end
    
    
    if strcmp(stimdata(sweepnum).Amplifiermode,'C-Clamp') & strcmp(rawdata(sweepnum).channellabel(1:end-2),'Vmon')
        stimstep=ceil(RSrisestep/2);
        stimdata(sweepnum).yforbridge=[stimdata(sweepnum).y(1:stimstep),stimdata(sweepnum).y(1:end-stimstep)];
        if RSrisestep==1
            stimdata(sweepnum).yforbridge=stimdata(sweepnum).yforbridge';
        else
            stimdata(sweepnum).yforbridge=moving(stimdata(sweepnum).yforbridge,RSrisestep);
        end
        eltolas=ceil(RSrisestep/2);
        stimdata(sweepnum).yforbridge=[ones(eltolas,1)*stimdata(sweepnum).yforbridge(1);stimdata(sweepnum).yforbridge(1:end-eltolas)];
        if size(stimdata(sweepnum).yforbridge,1)>size(stimdata(sweepnum).yforbridge,2)
            stimdata(sweepnum).yforbridge=stimdata(sweepnum).yforbridge';
        end
        RShs=find(abs(diff(stimdata(sweepnum).y))>0);
        RShs(find(RShs<10 |RShs>length(stimdata(sweepnum).y)-10))=[];
        Rallhs=[1,RShs,length(stimdata(sweepnum).y)];
        stimd=diff(stimdata(sweepnum).y);
        RSstimdiffs=stimd(RShs);
        RSbaseh1=[];
        RSbaseh2=[];
        for i=1:length(RShs)
            RSbaseh1(i,:)=RShs(i)-RSbaselinestep-1:RShs(i)-1;
            RSbaseh2(i,:)=RShs(i)+RSrisestep+1:RShs(i)+RSrisestep+RSbaselinestep+1;
            RS(i)=(mean(rawdata(sweepnum).y(RSbaseh2(i,:)))-mean(rawdata(sweepnum).y(RSbaseh1(i,:))))/RSstimdiffs(i);
            Rall(i)=(mean(rawdata(sweepnum).y(Rallhs(i+1):Rallhs(i+2)))-mean(rawdata(sweepnum).y(Rallhs(i):Rallhs(i+1))))/RSstimdiffs(i);
            RStime(i)=rawdata(sweepnum).realtime+time(RShs(i));
        end
        
        if any(Rall<10*10^6); %Hogyha RS+Rin kisebb, mint 10 MOhm, akkor valszeg ott nincs is stimulus..
%             figure(112)
%             clf
%             subplot(2,1,1)
%             plot(rawdata(sweepnum).y);
%             hold on;
%             plot(RSbaseh1,rawdata(sweepnum).y(RSbaseh1),'ro');
%             plot(RSbaseh2,rawdata(sweepnum).y(RSbaseh2),'go');
%             subplot(2,1,2)
%             plot(stimdata(sweepnum).y);
%             hold on;
%             plot(RSbaseh1,stimdata(sweepnum).y(RSbaseh1),'ro');
%             plot(RSbaseh2,stimdata(sweepnum).y(RSbaseh2),'go');
%             
%             
%             pause;
            stimdata(sweepnum).y=stimdata(sweepnum).y*0+rawdata(sweepnum).Amplifierholding(preamplnum);
            stimdata(sweepnum).yforbridge=stimdata(sweepnum).y;
        else
            %500pA-nél nagyobb áramokkal nem foglalkozunk
            rstokeep=find(abs(RSstimdiffs)<500*10^-12);
            RSs=[RSs,RS(rstokeep)];
            RStimes=[RStimes,RStime(rstokeep)];
            Ralls=[Ralls,Rall(rstokeep)];
        end
        if plotRSvalues==2
            
            figure(1)
            clf
            subplot(3,1,1)
            plot(time,rawdata(sweepnum).y)
            hold on
            plot(time(stimdata(sweepnum).segmenths),rawdata(sweepnum).y(stimdata(sweepnum).segmenths),'ro')
            plot(time(RSbaseh1),rawdata(sweepnum).y(RSbaseh1),'rx')
            plot(time(RSbaseh2),rawdata(sweepnum).y(RSbaseh2),'gx')
            subplot(3,1,2)
            plot(time,stimdata(sweepnum).y)
            ylim([-950*10^-12 950*10^-12])
            disp(RS)
            subplot(3,1,3)
            
            pause
        end
    else
        stimdata(sweepnum).yforbridge=stimdata(sweepnum).y*0;
    end
end
outlieridx=find(RSs>nanmean(RSs)+3*nanstd(RSs) | RSs<nanmean(RSs)-3*nanstd(RSs) | RSs<5*10^6 | RSs > 100*10^6);
RSs(outlieridx)=[];
RStimes(outlieridx)=[];
Ralls(outlieridx)=[];
outlieridx=find(RSs>nanmean(RSs)+3*nanstd(RSs) | RSs<nanmean(RSs)-3*nanstd(RSs) | RSs<0 | RSs > 100*10^6);
RSs(outlieridx)=[];
RStimes(outlieridx)=[];
Ralls(outlieridx)=[];
for sweepnum=1:length(rawdata)
    idx=RStimes<=rawdata(sweepnum).realtime+poolRStime &  RStimes>=rawdata(sweepnum).realtime-poolRStime;
    stimdata(sweepnum).RS=nanmean(RSs(idx));
end
if plotRSvalues==1
    figure(2323)
    clf
    subplot(2,1,1)
    plot(RStimes,RSs,'ko')
    hold on
    plot([rawdata.realtime],[stimdata.RS],'ro')
    subplot(2,1,2)
    plot(RStimes,Ralls,'ko')
end
for sweepnum=1:length(rawdata)
    if isnan(stimdata(sweepnum).RS)
        clear trs
        if sweepnum>1
            trs(1)=stimdata(sweepnum-1).RS;
        else
            trs(1)=NaN;
        end
        idx=find(~isnan([stimdata(sweepnum:end).RS]),1,'first');
        if isempty(idx)
            trs(2)=NaN;
        else
            trs(2)=stimdata(idx+sweepnum-1).RS;
        end
        stimdata(sweepnum).RS=nanmean(trs);
    end 
end
for sweepnum=1:length(rawdata) % if there is no RS measurement, RS is 0
    if isnan(stimdata(sweepnum).RS)
        stimdata(sweepnum).RS=0;
    end
end
if plotRSvalues==1
    subplot(2,1,1)    
    plot([rawdata.realtime],[stimdata.RS],'bx')
    xlabel('time (s)')
    ylabel('RS (Ohms)')
    pause
end

%creating PGF data, calculating RS