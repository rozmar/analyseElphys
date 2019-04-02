function aE_PSD_export(dirs,xlsdata,valtozok)
if ~isfield(valtozok,'onlyawake')
    valtozok.onlyawake = 0;
end
if isfield(valtozok,'breathing') & valtozok.breathing ==1
    sourcedir = dirs.breathingdir;
    destdir = dirs.breathingdir_psd;
    variablename = 'rawdata';
elseif isfield(valtozok,'high') & valtozok.high ==1
    sourcedir =dirs.bridgeddir;
    destdir = dirs.PSDdir_high;
    variablename = 'bridgeddata';
    valtozok.breathing=0;
elseif isfield(valtozok,'log') & valtozok.log ==1
    sourcedir =dirs.bridgeddir;
    destdir = dirs.PSDdir_log;
    variablename = 'bridgeddata';
    valtozok.breathing=0;
else
    sourcedir =dirs.bridgeddir;
    destdir = dirs.PSDdir;
    variablename = 'bridgeddata';
    valtozok.breathing=0;
end
parameters=valtozok.parameters;
for xlsi=length(xlsdata):-1:1
    if isfield(xlsdata,'anaesthesia')
        isitawake=any(strfind(xlsdata(xlsi).anaesthesia,'awake'));
    else
        isitawake=0;
    end
    if (valtozok.analyseonlyfield==0 | xlsdata(xlsi).field==1) & ( valtozok.onlyawake==0 | (valtozok.onlyawake==1 & isitawake==1))
        fname=[xlsdata(xlsi).HEKAfname,'_',xlsdata(xlsi).G_S_C];
        a=dir([destdir,fname,'.mat']);
        if isfield(valtozok,'breathing') & isfield(valtozok,'Thermosensor_channel')
            if valtozok.breathing ==1 & ~isnan(xlsdata(xlsi).Thermosensor_channel)
                cangoon=true;
            else
                cangoon=false;
            end
        else
            cangoon=true;
        end
        if (isempty(a) | valtozok.overwrite==1) & cangoon 
            temp=load([sourcedir,fname],variablename);
            bridgeddata=temp.(variablename);
            PSDdata=struct;
            for sweepi=1:length(bridgeddata)
                if length(bridgeddata(sweepi).y)*bridgeddata(sweepi).si>.5
                    
                    if ischar(valtozok.downsamplenum)
                        %%
                        si=bridgeddata(sweepi).si;
                        cutoffreq=parameters.max*10;
                        [b,a]=butter(1,cutoffreq/(1/si)/2,'low');
                        y = filtfilt(b,a,bridgeddata(sweepi).y);
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
%                     %%
%                     window=round(10/newsi); % 2 s sliding window
%                     noverlap=round(window-1); % downsampled step size
%                     [pm,fv,~] = spectrogram(y,window,noverlap,[parameters.min:parameters.step:parameters.max] ,1/newsi);
                    %%
                    PSDdata(sweepi).y=y;
                    PSDdata(sweepi).powerMatrix=powerMatrix;
                    PSDdata(sweepi).frequencyVector=frequencyVector;
                    PSDdata(sweepi).si_powerMatrix=newsi;
                    PSDdata(sweepi).realtime=bridgeddata(sweepi).realtime;
                    PSDdata(sweepi).timertime=bridgeddata(sweepi).timertime;
%                                     %%
%                                     parameters.wavenumber=9;
%                                     parameters.waveletlength=60;
%                                     parameters.taperlength=60;
%                                     [powerMatrix, frequencyVector]  = calculateTFDPower((y)', parameters);
%                                     figure(7)
%                                     clf
%                                     subplot(3,1,1)
%                                     plot([1:length(y)]*parameters.interval,y)
%                                     axis tight
%                                     subplot(3,1,2)
%                                     imagesc([1:length(y)]*parameters.interval+bridgeddata(sweepi).realtime,frequencyVector,powerMatrix)
%                                     set(gca,'YDir','normal');
%                                     caxis([0 .001]);
%                                     subplot(3,1,3)
%                                     imagesc([1:length(y)]*parameters.interval+bridgeddata(sweepi).realtime,fv,abs(pm))
%                                     set(gca,'YDir','normal');
%                                     caxis([0 .09]);
%                     %                 pause
%                     % %
                else
                    PSDdata(sweepi).powerMatrix=[];
                    PSDdata(sweepi).frequencyVector=[];
                    PSDdata(sweepi).si_powerMatrix=[];
                    PSDdata(sweepi).realtime=bridgeddata(sweepi).realtime;
                    PSDdata(sweepi).timertime=bridgeddata(sweepi).timertime;
                end
                
            end
            %% compressing PSD data
            for sweepi=1:length(PSDdata)
                if ~isempty(PSDdata(sweepi).powerMatrix)
                    maxval=max(PSDdata(sweepi).powerMatrix(:));
                    minval=min(PSDdata(sweepi).powerMatrix(:));
                    range=maxval-minval;
                    szorzo=(2^32-1)/range;
                    PSDdata(sweepi).powerMatrix=uint32((PSDdata(sweepi).powerMatrix-minval)*szorzo);
                    PSDdata(sweepi).compress_offset=minval;
                    PSDdata(sweepi).compress_multiplier=1/szorzo;
                else
                    PSDdata(sweepi).compress_offset=[];
                    PSDdata(sweepi).compress_multiplier=[];
                end
            end
            %%
            save([destdir,fname],'PSDdata','-v7.3');
            disp(['PSD export from ',fname, ' is done'])
        else
            disp(['PSD export from ',fname, ' was already done'])
%             if ~isempty(a)
%                 load([destdir,fname]);
%                 if ~isfield(PSDdata,'compress_offset')
%                     for sweepi=1:length(PSDdata)
%                         if ~isempty(PSDdata(sweepi).powerMatrix)
%                             maxval=max(PSDdata(sweepi).powerMatrix(:));
%                             minval=min(PSDdata(sweepi).powerMatrix(:));
%                             range=maxval-minval;
%                             szorzo=(2^32-1)/range;
%                             PSDdata(sweepi).powerMatrix=uint32((PSDdata(sweepi).powerMatrix-minval)*szorzo);
%                             PSDdata(sweepi).compress_offset=minval;
%                             PSDdata(sweepi).compress_multiplier=1/szorzo;
%                         else
%                             PSDdata(sweepi).compress_offset=[];
%                             PSDdata(sweepi).compress_multiplier=[];
%                         end
%                     end
%                     save([destdir,fname],'PSDdata','-v7.3');
%                     disp(['PSD export from ',fname, ' is compressed'])
%                 end
%             end
        end
        
        
    end
end