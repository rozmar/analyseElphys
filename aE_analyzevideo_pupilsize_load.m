function aE_analyzevideo_pupilsize_load(dirs,xlsdata)
savedir=[dirs.videodir,'eye/'];
for xlsi=1:length(xlsdata)
    setupname=xlsdata(xlsi).setup;
    filename=xlsdata(xlsi).HEKAfname(1:6);
    pupilsizedir=[locations.tgtardir,'VIDEOdata/',setupname,'/pupilsize/'];
    files=dir(pupilsizedir);
    files([files.isdir])=[];
    a=dir([savedir,xlsdata(xlsi).HEKAfname,'.mat']);
    if isempty(a)
        if length(files)>0
            for i=1:length(files)
                hyps=strfind(files(i).name,'_');
                files(i).filename=[files(i).name(hyps(1)-2:hyps(1)-1),files(i).name(hyps(1)+1:hyps(2)-1),files(i).name(hyps(2)+1:hyps(3)-1)];
                files(i).timestamp=[files(i).name(1:hyps(end)+2)];
                files(i).ext=files(i).name(end-2:end);
            end
            neededfiles=find(strcmp({files.filename},filename));
            if length(neededfiles)>0
                time=[];
                diameter=[];
                %             flash=[];
                for i=1:length(neededfiles)
                    %                 %%old method
                    %                 txt=load([pupilsizedir,files(neededfiles(i)).name]);
                    %                  time=[time;txt(:,3)];
                    %                 diameter=[diameter;txt(:,4)];
                    %                 flash=[flash;txt(:,5)];
                    %% new method - string read
                    fid=fopen([pupilsizedir,files(neededfiles(i)).name],'r');
                    txt=textscan(fid,'%s');
                    txt=txt{1};
                    fclose(fid);
                    dottedtimeidx = find(~cellfun(@isempty,regexp(txt,':..:')));
                    timenow=txt(dottedtimeidx+1);
                    diameternow=txt(dottedtimeidx+2);
                    timenow_num=zeros(size(timenow));
                    diameternow_num=zeros(size(diameternow));
                    for idxi=1:length(timenow)
                        if any(strfind(timenow{idxi},','))
                            timenow{idxi}(strfind(timenow{idxi},','))='.';
                        end
                        timenow_num(idxi)=str2num(timenow{idxi});
                        if any(strfind(diameternow{idxi},','))
                            diameternow{idxi}(strfind(diameternow{idxi},','))='.';
                        end
                        diameternow_num(idxi)=str2num(diameternow{idxi});
                    end
                    
                    time=[time;timenow_num];
                    diameter=[diameter;diameternow_num];
                    %%
                    
                end
                [time,ix]=sort(time);
                diameter=diameter(ix);
                %             flash=flash(ix);
                pupildata.diameter=diameter;
                pupildata.time=time;
                %             pupildata.flash=flash;

                save([savedir,xlsdata(xlsi).HEKAfname],'pupildata');
                disp([xlsdata(xlsi).HEKAfname,' - pupil diameter extraction is done'])
            end
        end
    end
end