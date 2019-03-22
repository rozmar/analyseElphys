function APwaves=persistent_generatefigures_generate_APwaves(dirs,xlsdata,xlsidx,valtozok)
% ID='1702021rm_6_1_3';
% xlsidx=find(strcmp({xlsdata.ID},ID));
% %%
% valtozok=struct;
% valtozok.debugmode=0;
% valtozok.xlimits=[0 60]+ 65104.1;
% valtozok.filter='gauss';% 'gauss' - gaussian filter, 'boxcar' - moving average
% valtozok.movingt=40;%microseconds for filtering: sigma in case of gaussian filter, width in case of boxcar
% valtozok.zerotime='threshold'; %threshold, apmaximum
% valtozok.timeback=.002; %ms back from zero time
% valtozok.timeforward=.002; %ms forward from zero time


ID=xlsdata(xlsidx).ID;
load([dirs.eventdir,'sorted/',ID,'.mat'],'eventdata');
load([dirs.bridgeddir,ID,'.mat'],'bridgeddata','stimdata');
neededevents=find(strcmpi({eventdata.type},'ap') & [eventdata.maxtime]>=valtozok.xlimits(1) & [eventdata.maxtime]<=valtozok.xlimits(2));

APwaves=struct;
apidxes=neededevents;
prevsweepnum=NaN;
for apii=1:length(apidxes)
    api=apidxes(apii);
    sweepnum=eventdata(api).sweepnum;
    if sweepnum~=prevsweepnum
        RS=stimdata(sweepnum).RS/10^6;
        y=bridgeddata(sweepnum).y;
        si=round(bridgeddata(sweepnum).si*10^6)/10^6;
        if valtozok.movingt>1
            if strcmp(valtozok.filter,'boxcar')
                valtozok.movingn=ceil(valtozok.movingt/si/1e6);
                yfiltered=moving(y,valtozok.movingn);
            elseif strcmp(valtozok.filter,'gauss')
                valtozok.movingn=ceil(valtozok.movingt/si/1e6)*10;
                hgauss=fspecial('gaussian',[valtozok.movingn,1],valtozok.movingt/si/1e6);
                yfiltered=imfilter(y,hgauss','replicate')';
            end    
        else
            yfiltered=y';
        end
        dyfiltered=diff(yfiltered)/si;
        if valtozok.movingt>1
            if strcmp(valtozok.filter,'boxcar')
                ddyfiltered=diff(moving(dyfiltered,valtozok.movingn));
                dddyfiltered=diff(moving(ddyfiltered,valtozok.movingn));
                ddddyfiltered=diff(moving(dddyfiltered,valtozok.movingn));
            elseif strcmp(valtozok.filter,'gauss')
                ddyfiltered=diff(imfilter(dyfiltered,hgauss,'replicate'));
                dddyfiltered=diff(imfilter(ddyfiltered,hgauss,'replicate'));
                ddddyfiltered=diff(imfilter(dddyfiltered,hgauss,'replicate'));
            end
        else
            ddyfiltered=diff(dyfiltered);
            dddyfiltered=diff(ddyfiltered);
            ddddyfiltered=diff(dddyfiltered);
        end
        dyfiltered=nanmean([[NaN;dyfiltered],[dyfiltered;NaN]],2);
        ddyfiltered=[ddyfiltered(1);ddyfiltered;ddyfiltered(end)];
        dddyfiltered=nanmean([[NaN;NaN;NaN;dddyfiltered],[dddyfiltered;NaN;NaN;NaN]],2);
        ddddyfiltered=[ddddyfiltered(1);ddddyfiltered(1);ddddyfiltered;ddddyfiltered(end);ddddyfiltered(end)];
        prevsweepnum=sweepnum;
        stepback=round(valtozok.timeback/si);
        stepforward=round(valtozok.timeforward/si);
    end
    
    if strcmp(valtozok.zerosettime,'threshold')
        centerh=eventdata(api).threshh;
    elseif strcmp(valtozok.zerosettime,'apmaximum')
        centerh=eventdata(api).maxh;
    end
    
    
    t=[-stepback:stepforward]*si*1000;

    if length(y)-stepforward<centerh
        vegh=length(y);
        nanvegen=stepforward-(length(y)-centerh);
    else
        vegh=centerh+stepforward;
        nanvegen=0;
    end
    if stepback>=centerh
        kezdeth=1;
        nanelejen=stepback-centerh+1;
    else
        kezdeth=centerh-stepback;
        nanelejen=0;
    end
    
    v=[nan(nanelejen,1);y(kezdeth:vegh)';nan(nanvegen,1)]';
    dv=[nan(nanelejen,1);dyfiltered(kezdeth:vegh);nan(nanvegen,1)]';%dyfiltered(threshh:maxh);
    ddv=[nan(nanelejen,1);ddyfiltered(kezdeth:vegh);nan(nanvegen,1)]';%ddyfiltered(threshh:maxh);
    dddv=[nan(nanelejen,1);dddyfiltered(kezdeth:vegh);nan(nanvegen,1)]';%dddyfiltered(threshh:maxh);
    ddddv=[nan(nanelejen,1);dddyfiltered(kezdeth:vegh);nan(nanvegen,1)]';%ddddyfiltered(threshh:maxh);

    APwaves(apii).t=t';
    APwaves(apii).v=v'*1000;
    APwaves(apii).dv=dv';
    APwaves(apii).ddv=ddv';
    APwaves(apii).dddv=dddv';
    APwaves(apii).ddddv=ddddv';
    APwaves(apii).si=si;
    APwaves(apii).RS=RS;
    APwaves(apii).maxtime=eventdata(api).maxtime;
    APwaves(apii).stimulated=eventdata(api).stimulated;
    APwaves(apii).axonalAP=eventdata(api).axonalAP;
    APwaves(apii).somaticAP=eventdata(api).somaticAP;
    APwaves(apii).threshv=eventdata(api).threshv;
    APwaves(apii).depolrate=eventdata(api).depolrate;
    APwaves(apii).baselineV=eventdata(api).baselineV;    
end