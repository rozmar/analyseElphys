function [sdvalue,sdvaluetime]=aE_calculatebaselineSD(time,y,valtozok,sdval)
Y=y;
TIME=time;
hossz=time(end);
segmentlength=3;
minsegmentlength=0;
segmentnum=ceil(hossz/segmentlength);


si=mode(diff(time));
step=round(valtozok.steptime/si);
stepback=ceil(step/2);
diffmovingstep=round(valtozok.diffmovingt/si);
filtmovingstep=round(valtozok.filtermovingtime/si);


for segmenti=1:(segmentnum)
%     if segmentnum>1
%         disp('upp')
%     end
    idx=(TIME>=(segmenti-1)*segmentlength & TIME<=segmenti*segmentlength);
    y=Y(idx);
    time=TIME(idx);
    if time(end)-time(1)>minsegmentlength
        %%
        sdvaluetime(segmenti)=time(round(length(time)/2));
        yfilt=imfilter(y,fspecial('average', [1,filtmovingstep]),'symmetric');
        ytemp=y;
        dy=diff(yfilt)/mode(diff(time));
        dyf=imfilter(dy,fspecial('average', [1,diffmovingstep]));
        tdyf=dyf;
        temptime=time(1:end-1);
        sdv=std(tdyf);
        %%
        while max(abs(tdyf))>sdval*sdv | max(abs(tdyf))>10 % törölgetjük a valószínűsíthető eseményeket, amíg van az adat deriváltjában egy adott érték feletti érték, vagy 5SD-n kívül van bármely derivált érték
            [~,loc]=max(abs(tdyf));
            if tdyf(loc)<0
                loc1=loc-1;
                if loc1<=0
                    loc1=1;
                end
                while loc1>1 & tdyf(loc1)<0
                    loc1=loc1-1;
                end
                loc2=loc+1;
                if loc2>=length(tdyf)
                    loc2=length(tdyf);
                end
                while loc2<length(tdyf) & tdyf(loc2)<0
                    loc2=loc2+1;
                end
            else
                loc1=loc-1;
                if loc1<=0
                    loc1=1;
                end
                while loc1>1 & tdyf(loc1)>0
                    loc1=loc1-1;
                end
                loc2=loc+1;
                if loc2>=length(tdyf)
                    loc2=length(tdyf);
                end
                while loc2<length(tdyf) & tdyf(loc2)>0
                    loc2=loc2+1;
                end
            end
            tdyf(loc1:loc2)=[];
            temptime(loc1:loc2)=[];
            sdv=std(tdyf);
        end
        sdvalue(segmenti)=sdv*valtozok.eventminsdval;
        
    end
end
%%
% dyf=dy;
% sdszorzo=3;
% sdwindow=.005;
% sdstep=round(sdwindow/si);
% dysd=moving(dyf,sdstep,'std');
% dymean=moving(dyf,sdstep,'mean');
% dysd=[nan(round(sdstep/2),1);dysd];
% dysd=dysd(1:length(dyf));
% dymean=[nan(round(sdstep/2),1);dymean];
% dymean=dymean(1:length(dyf));
% 
% eventidx=find(dyf'>dymean+dysd*sdszorzo | dyf'<dymean-dysd*sdszorzo);
% 
% figure(1)
% clf
% h(1)=subplot(2,1,1);
% hold on
% plot(dyf,'k-')
% plot(dymean,'k-','LineWidth',2)
% plot(dymean+dysd*sdszorzo,'r-')
% plot(dymean-dysd*sdszorzo,'b-')
% h(2)=subplot(2,1,2);
% plot(yfilt)
% hold on
% plot(eventidx,yfilt(eventidx),'ko')
% linkaxes(h,'x')
% %%
% 
% disp('a')