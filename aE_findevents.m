function aE_findevents(valtozok,dirs,parallelcount)
sdval=5;
files=dir(dirs.bridgeddir);
files([files.isdir])=[];
delete([dirs.eventparaleldir,'*.mat']);
for filenum=length(files):-1:1%1:length(files) % végigmegyünk az összes file-n
    a=dir([dirs.eventdir,files(filenum).name]);
    if isempty(a) || valtozok.overwrite==1
        workingmatlabnum=parallelcount+5;
        while workingmatlabnum>parallelcount-1
            filesunderprogress=dir(dirs.eventparaleldir);
            filesunderprogress([filesunderprogress.isdir])=[];
            workingmatlabnum=length(filesunderprogress);
        end
        source=[dirs.eventparaleldir,files(filenum).name];
        destination=[dirs.eventdir,files(filenum).name];
        save(source,'valtozok','dirs','files','filenum','sdval')
        pause(1);
%         aE_findevents_core(source,destination);
        unix(['matlab -nosplash -nodisplay -nojvm -r "aE_findevents_core(''',source,''',''',destination,'''), quit()"&']);
%         disp(['matlab -nosplash -nodisplay -nojvm -r "aE_findevents_core(''',source,''',''',destination,'''), quit()"&']);
        disp(['eventfinding started on',files(filenum).name]);
    else
        disp([files(filenum).name,' eventfinding already done.. skipped'])
    end
    
end
end