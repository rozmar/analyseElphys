function aE_cross_spectrum_export(dirs,xlsdata,valtozok)


sourcedir1 =dirs.bridgeddir;
label1='field';
sourcedir2 =dirs.breathingdir;
label2='breathing';
destdir = dirs.cross_spectrum_breathing;
variablename1 = 'bridgeddata';
variablename2 = 'rawdata';

parameters=valtozok.parameters;
for xlsi=length(xlsdata):-1:1
    if (valtozok.analyseonlyfield==0 | xlsdata(xlsi).field==1)
        fname=[xlsdata(xlsi).HEKAfname,'_',xlsdata(xlsi).G_S_C];
        a=dir([destdir,fname,'.mat']);
        if (isempty(a) | valtozok.overwrite==1) & ~isnan(xlsdata(xlsi).Thermosensor_channel)
            temp=load([sourcedir1,fname],variablename1);
            data1=temp.(variablename1);
            temp=load([sourcedir2,fname],variablename2);
            data2=temp.(variablename2);
            CROSSdata=struct;
            %%
            for sweep1num=1:length(data1)
                sweep2num=find(data1(sweep1num).realtime==[data2.realtime]);
                if ~isempty(sweep2num)
                    if isempty(fieldnames(CROSSdata));
                        NEXT=1;
                    else
                        NEXT=length(CROSSdata)+1;
                    end
                     if ischar(valtozok.downsamplenum)
                        si=data1(sweep1num).si;
                        cutoffreq=parameters.max*10;
                        [b,a]=butter(1,cutoffreq/(1/si)/2,'low');
                        y1 = filtfilt(b,a,data1(sweep1num).y);
                        y1=downsample(y1,round((1/cutoffreq)/si));
                        y2 = filtfilt(b,a,data2(sweep2num).y);
                        y2=downsample(y2,round((1/cutoffreq)/si));
                        newsi=1/cutoffreq;
                        parameters.interval=newsi;
                    else
                        parameters.interval=bridgeddata(sweepi).si*valtozok.downsamplenum;% - sample interval for the signal
                        y1=data1(sweep1num).y;
                        y1=downsample(y1,valtozok.downsamplenum);
                        y2=data2(sweep2num).y;
                        y2=downsample(y2,valtozok.downsamplenum);
                        newsi=data1(sweep1num).si*valtozok.downsamplenum;
                     end
                    
                    [powerMatrix1, frequencyVector1, CoeffMatrix1, ~]  = calculateTFDPower(y1', parameters);
                    [powerMatrix2, frequencyVector2, CoeffMatrix2, ~]  = calculateTFDPower(y2', parameters);
                    CROSSdata(NEXT).y1=y1;
                    CROSSdata(NEXT).y2=y2;
                    CROSSdata(NEXT).y1_label=label1;
                    CROSSdata(NEXT).y2_label=label2;
                    
                    %%
%                     h = fspecial('gaussian', [101,101],[30,30]);
                    h = hamming(round(4/newsi))';
                    Cross_spectrum=CoeffMatrix1.*conj(CoeffMatrix2);
%                     coherence=abs(Cross_spectrum).^2./(abs(CoeffMatrix1).^2.*abs(CoeffMatrix2).^2);
                    coherence_f=abs(imfilter(Cross_spectrum,h)).^2./imfilter(abs(CoeffMatrix1).^2,h).*imfilter(abs(CoeffMatrix2).^2,h);
                    phase=atan(imag(Cross_spectrum)./real(Cross_spectrum));
%                     %%
%                     figure(1)
%                     subplot(3,2,1)
%                     imagesc(powerMatrix1)
%                     subplot(3,2,2)
%                     imagesc(imfilter(powerMatrix1,h))
%                     subplot(3,2,3)
%                     imagesc(powerMatrix2)
%                     subplot(3,2,4)
%                     imagesc(imfilter(powerMatrix2,h))
%                     subplot(3,2,5)
%                     imagesc(coherence_f)
%                     caxis([0 .01])
%                     subplot(3,2,6)
%                     imagesc(phase)
% %                     pause
%                     caxis([0 .01])
                    %%
                    CROSSdata(NEXT).coherence=coherence_f;
                    CROSSdata(NEXT).phase=phase;
                    CROSSdata(NEXT).frequencyVector=frequencyVector1;
                    CROSSdata(NEXT).si_cross=newsi;
                    CROSSdata(NEXT).realtime=data1(sweep1num).realtime;
                    CROSSdata(NEXT).timertime=data1(sweep1num).timertime;
                end
                
%                 if length(bridgeddata(sweepi).y)*bridgeddata(sweepi).si>.5
%                     
%                     
% %                     %%
% %                     window=round(10/newsi); % 2 s sliding window
% %                     noverlap=round(window-1); % downsampled step size
% %                     [pm,fv,~] = spectrogram(y,window,noverlap,[parameters.min:parameters.step:parameters.max] ,1/newsi);
%                     %%
%                     PSDdata(sweepi).y=y;
%                     PSDdata(sweepi).powerMatrix=powerMatrix;
%                     PSDdata(sweepi).frequencyVector=frequencyVector;
%                     PSDdata(sweepi).si_powerMatrix=newsi;
%                     PSDdata(sweepi).realtime=bridgeddata(sweepi).realtime;
%                     PSDdata(sweepi).timertime=bridgeddata(sweepi).timertime;
% %                                     %%
% %                                     parameters.wavenumber=9;
% %                                     parameters.waveletlength=60;
% %                                     parameters.taperlength=60;
% %                                     [powerMatrix, frequencyVector]  = calculateTFDPower((y)', parameters);
% %                                     figure(7)
% %                                     clf
% %                                     subplot(3,1,1)
% %                                     plot([1:length(y)]*parameters.interval,y)
% %                                     axis tight
% %                                     subplot(3,1,2)
% %                                     imagesc([1:length(y)]*parameters.interval+bridgeddata(sweepi).realtime,frequencyVector,powerMatrix)
% %                                     set(gca,'YDir','normal');
% %                                     caxis([0 .001]);
% %                                     subplot(3,1,3)
% %                                     imagesc([1:length(y)]*parameters.interval+bridgeddata(sweepi).realtime,fv,abs(pm))
% %                                     set(gca,'YDir','normal');
% %                                     caxis([0 .09]);
% %                     %                 pause
% %                     % %
%                 else
%                     PSDdata(sweepi).powerMatrix=[];
%                     PSDdata(sweepi).frequencyVector=[];
%                     PSDdata(sweepi).si_powerMatrix=[];
%                     PSDdata(sweepi).realtime=bridgeddata(sweepi).realtime;
%                     PSDdata(sweepi).timertime=bridgeddata(sweepi).timertime;
%                 end
                
            end
            %% compressing PSD data
            for sweepi=1:length(CROSSdata)
                if ~isempty(CROSSdata(sweepi).coherence)
                    maxval=max(CROSSdata(sweepi).coherence(:));
                    minval=min(CROSSdata(sweepi).coherence(:));
                    range=maxval-minval;
                    szorzo=(2^32-1)/range;
                    CROSSdata(sweepi).coherence=uint32((CROSSdata(sweepi).coherence-minval)*szorzo);
                    CROSSdata(sweepi).compress_offset=minval;
                    CROSSdata(sweepi).compress_multiplier=1/szorzo;
                else
                    CROSSdata(sweepi).compress_offset=[];
                    CROSSdata(sweepi).compress_multiplier=[];
                end
            end
            %%
            save([destdir,fname],'CROSSdata','-v7.3');
            disp(['coherence export from ',fname, ' is done'])
        else
            disp(['coherence export from ',fname, ' was already done'])
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