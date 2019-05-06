function aE_analyzevideo_normalize_pupil_movement(dirs,xlsdata)
videopercentiles=struct;
alldata=struct;
files=dir([dirs.videodir,'eye/']);
files([files.isdir])=[];
for filei=1:length(files)
    
    load([dirs.videodir,'ROIs/',files(filei).name]);
    neededroi=find(strcmp({ROIdata.ROIname},'Eye'));
    stats=regionprops(ROIdata(neededroi).mask,'MajorAxisLength','MinorAxisLength');
    alldata(filei).MajorAxisLength=stats.MajorAxisLength;
    alldata(filei).MinorAxisLength=stats.MinorAxisLength;
    
    load([dirs.videodir,'eye/',files(filei).name]);
    load([dirs.videodir,'movement/',files(filei).name],'videodata');
    xlsidx=find(strcmp({xlsdata.HEKAfname},files(filei).name(1:end-4)),1,'first');
    
    moviename=[videodata(end).timestamp,'.avi'];
    setupname=xlsdata(xlsidx).setup;
    locations=marcicucca_locations;
    moviefiletoplay=[locations.tgtardir,'VIDEOdata/',setupname,'/',moviename];
    movieobj=VideoReader([moviefiletoplay]);
    szorzo=movieobj.Height/size(videodata(1).originalpic_all,1);
    %%
    
    alldata(filei).pupildiameter=pupildata.diameter'/alldata(filei).MinorAxisLength/szorzo*2;
    pupildiametersnow=sort(alldata(filei).pupildiameter);
    percentile95=pupildiametersnow(round(.95*length(pupildiametersnow)));
    if percentile95>1
        alldata(filei).pupildiameter=alldata(filei).pupildiameter/percentile95;
    end
    alldata(filei).pupildiameter_raw=pupildata.diameter';
%     alldata(filei).pupildiameter=pupildata.diameter';
    alldata(filei).pupildiameter=moving(alldata(filei).pupildiameter,2)';
    alldata(filei).pupildiameter_raw=moving(alldata(filei).pupildiameter_raw,2)';

    
end
sorteddiameter=sort([alldata.pupildiameter]);
pupilpercentiles=zeros(100,1);
for i=1:100
    pupilpercentiles(i)=sorteddiameter(round(length(sorteddiameter)*(i/100)));
end
videopercentiles.pupilpercentiles=pupilpercentiles;

allmovdata=struct;
files=dir([dirs.videodir,'movement/']);
files([files.isdir])=[];
for filei=1:length(files)
    load([dirs.videodir,'movement/',files(filei).name],'videodata');
    for fieldi=1:length(videodata);
        if isempty(fieldnames(allmovdata))
            NEXT=1;
        else
            NEXT=length(allmovdata)+1;
        end
        if length(videodata(fieldi).movement_all)>2
            allmovdata(NEXT).movementall=moving(videodata(fieldi).movement_all,2)';
        end
    end
end
sortedmovdata=sort([allmovdata.movementall]);
movementpercentiles=zeros(100,1);
for i=1:100
    movementpercentiles(i)=sortedmovdata(round(length(sortedmovdata)*(i/100)));
end
videopercentiles.movementpercentiles=movementpercentiles;
save([dirs.videodir,'percentiles'],'videopercentiles');