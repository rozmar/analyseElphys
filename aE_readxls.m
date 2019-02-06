function xlsdata=aE_readxls(path)
% %% Read XLS
[~,~,xlsdataraw] = xlsread(path);
drugnameidx=[];
drugtimeidx=[];
drugwashtimeidx=[];
drugconcidx=[];
numericrows={'age','Channel','startT','endT','locationX','locationy','locationz','slicethickness','Ca','Mg'};
for i=1:size(xlsdataraw,2)
    varname=char(xlsdataraw(1,i));
    if any(any(strcmpi(numericrows,varname)))
        itshouldbenumeric=true;
    else
        itshouldbenumeric=false;
    end
    for j=2:size(xlsdataraw,1)
        value=cell2mat(xlsdataraw(j,i));
        if itshouldbenumeric & ischar(value) %any(strfind(varname,'T'))
            xlsdata(j-1).(varname)=str2num(value);
        else
           xlsdata(j-1).(varname)=value;
           if ischar(value) & strcmp(value,'NaN')
               xlsdata(j-1).(varname)=NaN;
           end
        end
         
    end
    if any(strfind(varname,'DrugName'))
        drugnameidx=[drugnameidx,i];
    elseif any(strfind(varname,'DrugWashinTime'))
        drugtimeidx=[drugtimeidx,i];
    elseif any(strfind(varname,'DrugConcentration'))
        drugconcidx=[drugconcidx,i];
	elseif any(strfind(varname,'DrugWashoutTime'))
        drugwashtimeidx=[drugwashtimeidx,i];
    end
end
varname='ID';
for j=2:size(xlsdataraw,1)
    if isfield(xlsdata,'HEKAfname')
        xlsdata(j-1).(varname)=[xlsdata(j-1).HEKAfname,'_',xlsdata(j-1).G_S_C];
    else
        xlsdata(j-1).(varname)=[xlsdata(j-1).AXONfname,'_',xlsdata(j-1).Cellname];
    end
    xlsdata(j-1).drugnum=0;
    for i=1:length(drugnameidx)
        if ~any(strfind(cell2mat(xlsdataraw(j,drugnameidx(i))),'NaN'))
            xlsdata(j-1).drugnum=xlsdata(j-1).drugnum+1;
            if any(strfind(xlsdata(j-1).(['DrugName',num2str(i)]),'Calcium'))
                xlsdata(j-1).(['DrugName',num2str(i)])=['Calcium_',xlsdata(j-1).(['DrugConcentration',num2str(i)])];
            end
            xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugName=xlsdata(j-1).(['DrugName',num2str(i)]);%cell2mat(xlsdataraw(j,drugnameidx(xlsdata(j-1).drugnum)));
            xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugWashinTime=cell2mat(xlsdataraw(j,drugtimeidx(xlsdata(j-1).drugnum)));
            xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugWashoutTime=cell2mat(xlsdataraw(j,drugwashtimeidx(xlsdata(j-1).drugnum)));
            xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugConcentration=cell2mat(xlsdataraw(j,drugconcidx(xlsdata(j-1).drugnum)));
            if any(strfind(xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugConcentration,'mM'))
                szam=str2num(xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugConcentration(1:length(xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugConcentration)-2));
            elseif any(strfind(xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugConcentration,'microM'))
                szam=str2num(xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugConcentration(1:length(xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugConcentration)-6))/1000;
            elseif any(strfind(xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugConcentration,'nanoM'))
                szam=str2num(xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugConcentration(1:length(xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugConcentration)-6))/1000000;
            elseif any(strfind(xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugConcentration,'NaN'))
                szam=nan;
            else
                disp(['error - Unknown unit in ', xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugConcentration])
                pause
            end
            if ischar(xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugWashinTime)
                xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugWashinTime=str2num(xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugWashinTime);
            end
            if ischar(xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugWashoutTime)
                xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugWashoutTime=str2num(xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugWashoutTime);
            end
            % a drugdata-ban az időket korrigáljuk, hogyha az éjféli
            % váltásnál elrontja a HEKA
            if xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugWashinTime<xlsdata(j-1).startT & xlsdata(j-1).startT>22*3600 & xlsdata(j-1).endT<2*3600 
                xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugWashinTime=xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugWashinTime+24*3600;
            end
            if xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugWashoutTime<xlsdata(j-1).startT & xlsdata(j-1).startT>22*3600 & xlsdata(j-1).endT<2*3600 
                xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugWashoutTime=xlsdata(j-1).drugdata(xlsdata(j-1).drugnum).DrugWashoutTime+24*3600;
            end
        end
    end
end

clear i j xlsdataraw varname
end