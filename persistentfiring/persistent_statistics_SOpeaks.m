%% looking for slow oscillations on the field
percentilestosave=[10:10:100];
plotit=0;
oscillationpeaks=struct;
window=5;%s
windowstep=2.5;
freqrange=[.6 4];
for xlsidx=1:length(xlsdata)
    if true%xlsdata(xlsidx).field==1
        ID=xlsdata(xlsidx).ID;
        a=dir([dirs.PSDdir,ID,'.mat']);
        peakdata=struct;
        if ~isempty(a)
            load([dirs.PSDdir,ID])
            
            
            for sweep=1:length(PSDdata)
                sweephossz=length(PSDdata(sweep).y)*PSDdata(sweep).si_powerMatrix;
                if sweephossz>=window*1.5
                    time=[1:length(PSDdata(sweep).y)]*PSDdata(sweep).si_powerMatrix-PSDdata(sweep).si_powerMatrix;
                    if isfield(PSDdata,'compress_offset')
                        PSDdata(sweep).powerMatrix=(double(PSDdata(sweep).powerMatrix)-PSDdata(sweep).compress_offset)*PSDdata(sweep).compress_multiplier;
                    end
                    PSDmedian=zeros(size(PSDdata(sweep).powerMatrix,1),round(sweephossz/windowstep)-1);
                    PSDmin=zeros(size(PSDdata(sweep).powerMatrix,1),round(sweephossz/windowstep)-1);
                    PSDmediantime=1:size(PSDmedian,2)*windowstep;
                    %             peakdata=struct;
                    for wini=1:round(sweephossz/windowstep)-1
                        [~,ettol]=min(abs(time-((wini)*windowstep-window/2)));
                        [~,eddig]=min(abs(time-((wini)*windowstep+window/2)));
                        PSDmedian(:,wini)=nanmedian(PSDdata(sweep).powerMatrix(:,ettol:eddig),2);
                        PSDmin(:,wini)=nanmin(PSDdata(sweep).powerMatrix(:,ettol:eddig),[],2);
                        PSDmedian=PSDmin;
                        [pks,locs,w,p]=findpeaks(PSDmedian(:,wini));
                        freqs=PSDdata(sweep).frequencyVector(locs);
                        needed=freqs>freqrange(1) & freqs<=freqrange(2);
                        pks=pks(needed);
                        locs=locs(needed);
                        w=w(needed);
                        p=p(needed);
                        [p,idx]=sort(p,'descend');
                        pks=pks(idx);
                        locs=locs(idx);
                        w=w(idx);
                        peaknum=2;
                        if length(locs)>peaknum+1
                            
                            idxnow=locs(1:peaknum);
                            if isempty(fieldnames(peakdata))
                                next=1;
                            else
                                next=length(peakdata)+1;
                            end
                            peakdata(next).peakval=pks(1:peaknum)';
                            peakdata(next).peakvalratio=pks(1:peaknum)'/pks(peaknum+1);
                            peakdata(next).peakwidth=w(1:peaknum)';
                            peakdata(next).prominence=p(1:peaknum)';
                            peakdata(next).prominenceratio=p(1:peaknum)'/p(peaknum+1);
                            peakdata(next).freq=PSDdata(sweep).frequencyVector(locs(1:peaknum));
                            peakdata(next).sweepnum=ones(size(peakdata(next).peakval))*sweep;
                        end
                        %                 figure(2)
                        %                 clf
                        %                 plot(PSDmedian(:,wini));
                        %                 hold on
                        %                 plot(idxnow,PSDmedian(idxnow,wini),'ro')
                        %                 pause
                    end
                    %             if any([peakdata.peakvalratio]>2)  %any([peakdata.peakval]>.001)
                    %                 %%
                    %                 figure(1)
                    %                 clf
                    %                 subplot(3,2,1)
                    %                 plot(time,PSDdata(sweep).y);
                    %                 subplot(3,2,2)
                    %                 plot(PSDdata(sweep).frequencyVector,PSDmedian);
                    %                 subplot(3,2,3)
                    %                 imagesc(time,PSDdata(sweep).frequencyVector,PSDdata(sweep).powerMatrix);
                    %                 set(gca,'YDir','normal');
                    %                 colormap linspecer
                    %                 ylabel('Frequency (Hz)')
                    %                 xlabel('Time (s)')
                    %                 subplot(3,2,4)
                    %                 imagesc(PSDmediantime,PSDdata(sweep).frequencyVector,PSDmedian);
                    %                 set(gca,'YDir','normal');
                    %                 colormap linspecer
                    %                 ylabel('Frequency (Hz)')
                    %                 xlabel('Time (s)')
                    %
                    %                 subplot(3,2,5)
                    %                 semilogy([peakdata.freq],[peakdata.peakval],'ko')
                    %                 xlabel('freq')
                    %                 ylabel('peak power')
                    %                 subplot(3,2,6)
                    %                 plot([peakdata.peakvalratio],[peakdata.freq],'ko')
                    %                 xlabel('peakvalratio')
                    %                 ylabel('freq')
                    %                 pause
                    %             end
                end
            end
        end
        if ~isempty(fieldnames(peakdata))
            if  any([peakdata.peakvalratio]>2) & plotit==1  %any([peakdata.peakval]>.001)
                %%
                
                
                [~,idx]=sort([peakdata.peakvalratio],'descend');
                
                peakvals=[peakdata.peakval];
                peakvals=peakvals(idx);
                
                sweepnums=[peakdata.sweepnum];
                sweepnums=sweepnums(idx);
                
                [sweepnums,ic,~]=unique(sweepnums,'stable');
                peakvals=peakvals(ic);
                
                sweep=sweepnums(1);
                peakval=peakvals(1);
                time=[1:length(PSDdata(sweep).y)]*PSDdata(sweep).si_powerMatrix;
                figure(1)
                clf
                subplot(3,2,1)
                plot(time,PSDdata(sweep).y);
                subplot(3,2,2)
                imagesc(time,PSDdata(sweep).frequencyVector,PSDdata(sweep).powerMatrix);
                set(gca,'YDir','normal');
                colormap linspecer
                caxis([0 peakval]);
                ylabel('Frequency (Hz)')
                xlabel('Time (s)')
                subplot(3,2,3)
                [n,xbin]=hist([peakdata.peakvalratio],[0:50]);
                bar(xbin,n/sum(n)*100)
                xlim([0 51])
                xlabel('peakvalratio')
                ylabel('percentage in bin')
                subplot(3,2,4)
                plot(xbin,(1-cumsum(n/sum(n)))*100)
                xlabel('peakvalratio')
                ylabel('percentage above value')
                subplot(3,2,5)
                semilogy([peakdata.freq],[peakdata.peakval],'ko')
                xlabel('freq')
                ylabel('peak power')
                subplot(3,2,6)
                plot([peakdata.peakvalratio],[peakdata.freq],'ko')
                xlabel('peakvalratio')
                ylabel('freq')
                
                figure(2)
                clf
                ennyisweep=min(length(sweepnums)-1,3);
                for i=1:ennyisweep
                    sweep=sweepnums(i+1);
                    peakval=peakvals(i+1);
                    time=[1:length(PSDdata(sweep).y)]*PSDdata(sweep).si_powerMatrix;
                    subplot(ennyisweep,2,(i-1)*2+1)
                    plot(time,PSDdata(sweep).y);
                    subplot(ennyisweep,2,(i-1)*2+2)
                    imagesc(time,PSDdata(sweep).frequencyVector,PSDdata(sweep).powerMatrix);
                    set(gca,'YDir','normal');
                    colormap linspecer
                    caxis([0 peakval]);
                    ylabel('Frequency (Hz)')
                    xlabel('Time (s)')
                end
                pause
                
            end
            %%
            [peakvalratios,idx]=sort([peakdata.peakvalratio],'ascend');
            
            peakratiopercentiles=zeros(size(percentilestosave));
            for perci=1:length(percentilestosave)
                idxnow=round(length(peakvalratios)*percentilestosave(perci)/100);
                idxnow=max(idxnow,1);
                idxnow=min(idxnow,length(peakvalratios));
                peakratiopercentiles(perci)=peakvalratios(idxnow);
            end
            
            peakvals=[peakdata.peakval];
            peakvals=peakvals(idx);
            
            if isempty(fieldnames(oscillationpeaks))
                NEXT=1;
            else
                NEXT=length(oscillationpeaks)+1;
            end
            oscillationpeaks(NEXT).ID=xlsdata(xlsidx).ID;
            for perci=1:length(percentilestosave)
                oscillationpeaks(NEXT).(['peakratio_',num2str(percentilestosave(perci)),'_percentile'])=peakratiopercentiles(perci);
            end
            
        else
            if isempty(fieldnames(oscillationpeaks))
                NEXT=1;
            else
                NEXT=length(oscillationpeaks)+1;
            end
            oscillationpeaks(NEXT).ID=xlsdata(xlsidx).ID;
            for perci=1:length(percentilestosave)
                oscillationpeaks(NEXT).(['peakratio_',num2str(percentilestosave(perci)),'_percentile'])=NaN;
                
            end
        end
        disp([ID,': slow oscillation peaks saved'])
    end
    
end
fieldek=fieldnames(oscillationpeaks);
for xlsnum=1:length(xlsdata) %% add field values to recordings
    fieldxlsnum=find(strcmp(xlsdata(xlsnum).HEKAfname,{xlsdata.HEKAfname}) & [xlsdata.field]==1);
    
    if isempty(fieldxlsnum)
        for fieldi=1:length(fieldek)
            fieldname=fieldek{fieldi};
            oscillationpeaks(xlsnum).([fieldname,'_field'])=NaN;
        end
    else
        fieldxlsnum=fieldxlsnum(1);
        for fieldi=1:length(fieldek)
            fieldname=fieldek{fieldi};
            oscillationpeaks(xlsnum).([fieldname,'_field'])=oscillationpeaks(fieldxlsnum).(fieldname);
        end
    end
end
save([dirs.basedir,'statistics/SOpeaks.mat'],'oscillationpeaks');