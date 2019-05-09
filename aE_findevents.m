function aE_findevents(valtozok,dirs,parallelcount,xlsdata)
% aE_findevents.m runs aE_findevents_core.m in a parallelized fashion.
% In the case of each processed file, a temporary .mat file is written in 
% dirs.eventparaleldir which contains the name of the file to be processed 
% and the valtozok variable. Then a new matlab instance is started to run 
% aE_findevents_core.m script for that given file. When the 
% aE_findevents_core.m finishes, it deletes the temporary .mat file in 
% dirs.eventparaleldir so another matlab instance can be started.
%
% The following parameters must be given: valtozok, dirs, parallelcount,
% xlsdata
%   -valtozok: this is the parameter set for aE_findevents_core.m. Only two
%   fields are relevant for this particular script:
%        -valtozok.overwrite: 1 if the already exported data should be
%        overwritten
%        -valtozok.overwritebefore: double in datenum format (datenum(datetime('today')))
%        if valtozok.overwrite==1, only files generated before this date
%        will be overwritten
%   -dirs: directory structure generated in analysElphys_main
%   -parallelcount: int - number of parallel matlab instances. 
%   -xlsdata: metadata of experiments generated in analysElphys_main
%
% See also aE_findevents_core, aE_calculatebaselineSD

sdval=valtozok.sdval;
files=dir(dirs.bridgeddir);
files([files.isdir])=[];
if parallelcount>1
     delete([dirs.eventparaleldir,'*.mat']);
end
for filenum=length(files):-1:1%filenum=1:length(files)%% végigmegyünk az összes file-n % 
    a=dir([dirs.eventdir,files(filenum).name]);
    if isfield(xlsdata,'ID')
        xlsnum=find(strcmp(files(filenum).name(1:end-4),{xlsdata.ID}));
    else
        xlsnum=NaN;
    end
    %%
    mehet=true;
    if isnan(xlsnum)
        mehet=false;
    end
    fieldstocheck={'juxta','field','axonal'};
    for fieldi=1:length(fieldstocheck)
        if ~isfield(xlsdata,fieldstocheck{fieldi});
            mehet=false;
        else
            if xlsdata(xlsnum).(fieldstocheck{fieldi})==1
                mehet=false;
            end
        end
    end
    %%
    if (isempty(a) || (valtozok.overwrite==1 & valtozok.overwritebefore>a.datenum)) & mehet %(isnan(xlsnum) | ((~isfield(xlsdata,'juxta ') | xlsdata(xlsnum).juxta==0) & (~isfield(xlsdata,'field') | xlsdata(xlsnum).field==0) & (~isfield(xlsdata,'axonal ') | xlsdata(xlsnum).axonal==0)))
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