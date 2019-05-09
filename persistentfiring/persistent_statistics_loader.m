function xlsdata=persistent_statistics_loader(dirs,xlsdata,projectdata)
statfilenames={'APstats.mat','SOpeaks.mat'};
statvarrnames={'APstatdata','oscillationpeaks'};
if projectdata.loadstatistics==1
    disp('loading statistics...')
    for stati=1:length(statfilenames)
        a=dir([dirs.basedir,'statistics/',statfilenames{stati}]);
        if ~isempty(a)
            temp=load([dirs.basedir,'statistics/',statfilenames{stati}]);
            statdata=temp.(statvarrnames{stati});
            
            xlsfields=fieldnames(xlsdata);
            fieldstoadd=fieldnames(statdata);
            todel=false(size(fieldstoadd));
            for i=1:length(fieldstoadd)
                if any(strcmp(fieldstoadd{i},xlsfields))
                    todel(i)=true;
                end
            end
            fieldstoadd(todel)=[];
            for xlsi=1:length(xlsdata)
                idx=find(strcmp(xlsdata(xlsi).ID,{statdata.ID}));
                if isempty(idx)
                    for fieldi=1:length(fieldstoadd)
                        xlsdata(xlsi).(fieldstoadd{fieldi})=NaN;
                    end
                else
                    for fieldi=1:length(fieldstoadd)
                        xlsdata(xlsi).(fieldstoadd{fieldi})=statdata(idx).(fieldstoadd{fieldi});
                    end
                end
                
            end
        end
    end
    
    fieldek=fieldnames(xlsdata);
for fieldi=1:length(fieldek)
    fieldnev=fieldek{fieldi};
    if any(cellfun(@isempty,{xlsdata.(fieldnev)}))%length({xlsdata.(fieldnev)})<length(xlsdata)
        [xlsdata((cellfun(@isempty,{xlsdata.(fieldnev)}))).(fieldnev)]=deal(NaN);
    end
end
end