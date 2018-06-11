function aE_PSD_export(dirs,xlsdata,valtozok);
parameters=valtozok.parameters;
for xlsi=length(xlsdata):-1:1
    
    if valtozok.analyseonlyfield==0 | xlsdata(xlsi).field==1
        fname=[xlsdata(xlsi).HEKAfname,'_',xlsdata(xlsi).G_S_C];
        a=dir([dirs.PSDdir,fname,'.mat']);
        if isempty(a)|valtozok.overwrite==1;
            load([dirs.bridgeddir,fname],'bridgeddata');
            PSDdata=struct;
            for sweepi=1:length(bridgeddata)
                if length(bridgeddata(sweepi).y)*bridgeddata(sweepi).si>.5
                    if ischar(valtozok.downsamplenum)
                        %%
                        si=bridgeddata(sweepi).si;
                        cutoffreq=parameters.max*10;
                        [b,a]=butter(1,cutoffreq/(1/si)/2,'low');
                        y = filter(b,a,bridgeddata(sweepi).y);
                        y=downsample(y,round((1/cutoffreq)/si));
                        newsi=1/cutoffreq;
                        parameters.interval=newsi;
                    else
                        parameters.interval=bridgeddata(sweepi).si*valtozok.downsamplenum;% - sample interval for the signal
                        y=bridgeddata(sweepi).y;
                        y=downsample(y,valtozok.downsamplenum);
                        newsi=bridgeddata(sweepi).si*valtozok.downsamplenum;
                    end
                    
                    [powerMatrix, frequencyVector]  = calculateTFDPower(y', parameters);
                    PSDdata(sweepi).y=y;
                    PSDdata(sweepi).powerMatrix=powerMatrix;
                    PSDdata(sweepi).frequencyVector=frequencyVector;
                    PSDdata(sweepi).si_powerMatrix=newsi;
                    PSDdata(sweepi).realtime=bridgeddata(sweepi).realtime;
                    PSDdata(sweepi).timertime=bridgeddata(sweepi).timertime;
                    %                 %%
                    %                 figure(1)
                    %                 clf
                    %                 subplot(2,1,1)
                    %                 plot(bridgeddata(sweepi).y)
                    %                 subplot(2,1,2)
                    %                 imagesc([1:length(y)]*parameters.interval+bridgeddata(sweepi).realtime,frequencyVector,powerMatrix)
                    %                 set(gca,'YDir','normal');
                    %
                    % %                 caxis([0 4]);
                    %                 pause
                    % %
                else
                    PSDdata(sweepi).powerMatrix=[];
                    PSDdata(sweepi).frequencyVector=[];
                    PSDdata(sweepi).si_powerMatrix=[];
                    PSDdata(sweepi).realtime=bridgeddata(sweepi).realtime;
                    PSDdata(sweepi).timertime=bridgeddata(sweepi).timertime;
                end
                
            end
            save([dirs.PSDdir,fname],'PSDdata','-v7.3');
            disp(['PSD export from ',fname, ' is done'])
        else
            disp(['PSD export from ',fname, ' was already done'])
        end
        
        
    end
end