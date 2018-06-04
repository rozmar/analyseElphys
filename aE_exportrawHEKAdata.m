function analysElphys_exportrawHEKAdata(dirs,xlsdata,overwrite)
for xlsidx=1:length(xlsdata)
    a=dir([dirs.rawexporteddir,xlsdata(xlsidx).ID,'.mat']);
    if isempty(a) | overwrite==1
        rawdata=HEKAexportbytime_main(xlsdata(xlsidx).HEKAfname,xlsdata(xlsidx).setup,xlsdata(xlsidx).Channel,xlsdata(xlsidx).startT,xlsdata(xlsidx).endT);
        if isfield(xlsdata,'Stimulation') & ~strcmp(xlsdata(xlsidx).Stimulation,'none') & ~strcmp(xlsdata(xlsidx).Stim_channel,'none')
            rawstimdata=HEKAexportbytime_main(xlsdata(xlsidx).HEKAfname,xlsdata(xlsidx).setup,xlsdata(xlsidx).Stim_channel,xlsdata(xlsidx).startT,xlsdata(xlsidx).endT);
            %%
            for sweepi=1:length(rawdata)
                if any(rawdata(sweepi).realtime==[rawstimdata.realtime])
                    idx=find(rawdata(sweepi).realtime==[rawstimdata.realtime]);
                    if any(strfind(xlsdata(xlsidx).Stimulation,'light'))
                        threshval=(max(rawstimdata(idx).y)-min(rawstimdata(idx).y))/2+min(rawstimdata(idx).y);
                        stimulation=rawstimdata(idx).y>threshval;
                    else
                        stimulation=rawstimdata(idx).y;
                    end
                else
                    stimulation=false(size(rawdata(sweepi).y));
                end
                rawdata(sweepi).stimulation=stimulation;
            end
        end
        save([dirs.rawexporteddir,xlsdata(xlsidx).ID],'rawdata','xlsdata','xlsidx','-v7.3')
        disp([xlsdata(xlsidx).ID,' done'])
    else
        disp([xlsdata(xlsidx).ID,' already done.. skipped'])
    end
end
xlsdataold=xlsdata;
end