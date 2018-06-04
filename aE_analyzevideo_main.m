%%
function aE_analyzevideo_main(valtozok,dirs)

locations=marcicucca_locations;
sampleinterval=valtozok.sampleinterval;
setupname=valtozok.setupname;
filename=valtozok.filename;

a=dir([dirs.videodir,'movement/',filename,'.mat']);
if isempty(a)|valtozok.overwrite==1
    
    rawvideodir=[locations.tgtardir,'VIDEOdata/',setupname,'/'];
%     pupilsizedir=[locations.tgtardir,'VIDEOdata/',setupname,'/pupilsize/'];
    files=dir(rawvideodir);
    files([files.isdir])=[];
    
    for i=1:length(files)
        hyps=strfind(files(i).name,'_');
        files(i).filename=[files(i).name(hyps(1)-2:hyps(1)-1),files(i).name(hyps(1)+1:hyps(2)-1),files(i).name(hyps(2)+1:hyps(3)-1)];
        files(i).timestamp=[files(i).name(1:hyps(end)+2)];
        files(i).ext=files(i).name(end-2:end);
    end
    neededvideos=find(strcmp({files.filename},filename(1:6))&strcmp({files.ext},'avi'));
    realfilei=0;
    for filei=1:length(neededvideos)
        disp(['extracting video info from file: ', files(neededvideos(filei)).name]);
        movieobj=VideoReader([rawvideodir,files(neededvideos(filei)).name]);
        correspondingtxt=find(strcmp(files(neededvideos(filei)).timestamp,{files.timestamp})&strcmp({files.ext},'txt'));
        nFrames = movieobj.NumberOfFrames;
        if  nFrames>10 & ~any(strfind(files(neededvideos(filei)).name,'eye'))
            realfilei=realfilei+1;
        %     txt=load([rawvideodir,files(correspondingtxt).name]);
        %     frametimes=txt(:,3);
%         %%
%         fid   = fopen([rawvideodir,files(correspondingtxt).name]);
%         clear texty
%         texty = textscan(fid,'%s%s%s%s%s');
%         fclose(fid);
        
        %%
        fid   = fopen([rawvideodir,files(correspondingtxt).name]);
        clear texty
        texty = textscan(fid,'%s');
        fclose(fid);
        %%
%         ezaz=0;
%         for i=1:length(texty)
%             a=texty{i}{1};
%             if any(regexp(a,':..:'))
%                 ezaz=i;
%             end
%         end
%         texty=texty{ezaz};
         a=cellfun(@any, regexp(texty{1},':..:'));
        texty=texty{1}(a);
            %%
        frametimes=zeros(size(texty));
        for i=1:length(texty)
            numnow=texty{i};
            points=strfind(numnow,':');
            if length(points)==2
                if any(strfind(numnow,','))
                    numnow(strfind(numnow,','))='.';
                end
                frametimes(i)=str2num(numnow(1:points(1)-1))*3600+str2num(numnow(points(1)+1:points(2)-1))*60+str2num(numnow(points(2)+1:end));
            else
                frametimes(i)=NaN;
            end
        end
        %%
        neededframes=1;
        while frametimes(neededframes(end))+sampleinterval<nanmax(frametimes) & neededframes(end)<length(frametimes)
            [~,nextidx]=min(abs(frametimes(neededframes(end)+1:end)-frametimes(neededframes(end))-sampleinterval));
            neededframes=[neededframes,neededframes(end)+nextidx];
        end
        neededframes(end)=[];
        
        %%
        nFrames = movieobj.NumberOfFrames;
        dummypic=rgb2gray(read(movieobj, 1));
        dummypic=imresize(dummypic,.1);
        vidHeight=size(dummypic,1);
        vidWidth=size(dummypic,2);
        clear mov
        mov=uint8(zeros(vidHeight,vidWidth,length(neededframes)));
        for k = 1:length(neededframes)
            mov(:,:,k) = imresize(rgb2gray(read(movieobj, neededframes(k))),.1);
            progressbar(neededframes(k)/nFrames);
        end
        %%
        close all
        neededframetimes=frametimes(neededframes);
        %     downmov=mov;%imresize(mov,ceil([size(mov,1)/10,size(mov,2)/10]));
        downmov=abs(diff(mov,1,3));
        neededframetimes=mean([neededframetimes(1:end-1),neededframetimes(2:end)],2);
        lineardownmov=reshape(downmov,size(downmov,1)*size(downmov,2),size(downmov,3));
        
        button=106;
        while button~=107
            figure(3)
            clf
            subplot(2,3,1)
            imagesc(mean(mov,3))
            subplot(2,3,2)
            bigpic=max(downmov,[],3);
            imagesc(bigpic);
            movingaveragebig=zeros(size(downmov,3),1);
            for framei=1:size(downmov,3)
                imnow=downmov(:,:,framei);
                movingaveragebig(framei)=mean(double(imnow(:)));
            end
            subplot(2,3,3)
            plot(neededframetimes,movingaveragebig);
            
             subplot(2,3,2)
%             [x,y]=ginput(2);
            x=[0 inf];
            y=[0 inf];
            xidxs=round([max(min(y),1):min(max(y),size(downmov,1))]);
            yidxs=round([max(min(x),1):min(max(x),size(downmov,2))]);
            
             subplot(2,3,4)
            imagesc(mean(mov(xidxs,yidxs,:),3))
            
            selectedpic=max(downmov(xidxs,yidxs,:),[],3);
            subplot(2,3,5)
            imagesc(selectedpic)
            movingaverageselected=zeros(size(downmov,3),1);
            for framei=1:size(downmov,3)
                imnow=downmov(xidxs,yidxs,framei);
                movingaverageselected(framei)=mean(double(imnow(:)));
            end
            subplot(2,3,6)
            plot(neededframetimes,movingaverageselected);
            
%             [~,~,button]=ginput(1);
            button=107;
        end
        
        
        %%
        %     [coeff,score]=princomp(double(lineardownmov'));
        %     %%
        %     close all
        %     princompstoshow=6;
        %     figure(10)
        %     clf
        %     for i=1:princompstoshow
        %         subplot(princompstoshow,2,(i-1)*2+1)
        %         plot((score(:,i)))
        %         subplot(princompstoshow,2,(i-1)*2+2)
        %         scorelayout(:,:,i)=reshape(coeff(:,i),size(downmov,1),size(downmov,2));
        %         imagesc(scorelayout(:,:,i))
        %         caxis([-.04 .04])
        %     end
        %%
        
        %     videodata(filei).scorelayout=scorelayout;
        %     videodata(filei).princomps=score(:,1:princompstoshow);
        videodata(realfilei).movement_all=movingaveragebig;
        videodata(realfilei).movement_selected=movingaverageselected;
        videodata(realfilei).originalpic_all=mean(mov,3);
        videodata(realfilei).diffpic_all=bigpic;
        videodata(realfilei).originalpic_selected=mean(mov(xidxs,yidxs,:),3);
        videodata(realfilei).diffpic_all=selectedpic;
        videodata(realfilei).time=neededframetimes;
        videodata(realfilei).si=sampleinterval;
        videodata(realfilei).timestamp=files(neededvideos(filei)).timestamp;
        rawvideo(realfilei).vid=mov;
        end
    end
    if exist('videodata','var')
        save([dirs.videodir,'movement/',filename],'videodata','rawvideo');
    end
else
    disp(['video analysis of ', filename,' is already done'])
end