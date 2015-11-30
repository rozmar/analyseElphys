function analysElphys_exportrawHEKAdata(dirs,xlsdata,overwrite)
for xlsidx=1:length(xlsdata)
    a=dir([dirs.rawexporteddir,xlsdata(xlsidx).ID,'.mat']);
    if isempty(a) | overwrite==1
        rawdata=HEKAexportbytime_main(xlsdata(xlsidx).HEKAfname,xlsdata(xlsidx).setup,xlsdata(xlsidx).Channel,xlsdata(xlsidx).startT,xlsdata(xlsidx).endT);
        save([dirs.rawexporteddir,xlsdata(xlsidx).ID],'rawdata','xlsdata','xlsidx')
        disp([xlsdata(xlsidx).ID,' done'])
    else
        disp([xlsdata(xlsidx).ID,' already done.. skipped'])
    end
end
xlsdataold=xlsdata;
end