function sdvalue=aE_calculatebaselineSD(time,y,valtozok)
si=mode(diff(time));
step=round(valtozok.steptime/si);
stepback=ceil(step/2);
diffmovingstep=round(valtozok.diffmovingt/si);
filtmovingstep=round(valtozok.filtermovingtime/si);

yfilt=imfilter(y,fspecial('average', [1,filtmovingstep]),'symmetric');
ytemp=y;
dy=diff(yfilt)/mode(diff(time));
dyf=imfilter(dy,fspecial('average', [1,diffmovingstep]));
tdyf=dyf;
temptime=time(1:end-1);
sdv=std(tdyf);
while max(abs(tdyf))>3*sdv | max(abs(tdyf))>10 % törölgetjük a valószínűsíthető eseményeket, amíg van az adat deriváltjában egy adott érték feletti érték, vagy 3SD-n kívül van bármely derivált érték
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
sdvalue=sdv*valtozok.eventminsdval;