function aE_findevents(valtozok,dirs,parallelcount,xlsdata)
sdval=3;
files=dir(dirs.bridgeddir);
files([files.isdir])=[];
if parallelcount>1
    delete([dirs.eventparaleldir,'*.mat']);
end
for filenum=1:length(files)% length(files):-1:1%% végigmegyünk az összes file-n %
    a=dir([dirs.eventdir,files(filenum).name]);
    if isfield(xlsdata,'ID')
        xlsnum=find(strcmp(files(filenum).name(1:end-4),{xlsdata.ID}));
    else
        xlsnum=NaN;
    end
    if (isempty(a) || (valtozok.overwrite==1 & valtozok.overwritebefore>a.datenum)) & (isnan(xlsnum)| ((~isfield(xlsdata,'juxta ') | xlsdata(xlsnum).juxta==0) & (~isfield(xlsdata,'field') | xlsdata(xlsnum).field==0) & (~isfield(xlsdata,'axonal ') | xlsdata(xlsnum).axonal==0)))
        %         workingmatlabnum=parallelcount+5;
        if parallelcount>1
            filesunderprogress=dir(dirs.eventparaleldir);
            filesunderprogress([filesunderprogress.isdir])=[];
            workingmatlabnum=length(filesunderprogress);
            while workingmatlabnum>parallelcount-1
                filesunderprogress=dir(dirs.eventparaleldir);
                filesunderprogress([filesunderprogress.isdir])=[];
                workingmatlabnum=length(filesunderprogress);
                pause(10);
            end
        end
        source=[dirs.eventparaleldir,files(filenum).name];
        destination=[dirs.eventdir,files(filenum).name];
        save(source,'valtozok','dirs','files','filenum','sdval')
        disp(['eventfinding started on ',files(filenum).name]);
        if parallelcount<=1
            aE_findevents_core(source,destination);
        else
            unix(['matlab -nosplash -nodisplay -nojvm -r "aE_findevents_core(''',source,''',''',destination,'''), quit()"&']);
        end
%         disp(['matlab -nosplash -nodisplay -nojvm -r "aE_findevents_core(''',source,''',''',destination,'''), quit()"&']);
        
    else
        disp([files(filenum).name,' eventfinding already done.. skipped'])
    end
    
end
end