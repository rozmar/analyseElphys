function aE_analyzevideo_PCA_on_ROIs(dirs,xlsdata,overwrite)
% This is an unfinished script. It works but the results are not
% satisfactory.


ROIdir=[dirs.videodir,'ROIs/'];
movementdir=[dirs.videodir,'movement/'];
files=dir(movementdir);
files([files.isdir])=[];
for filei=1:length(files)
    a=dir([ROIdir,files(filei).name]);
    if ~isempty(a)
    load([ROIdir,files(filei).name]);
    %     if true%~isfield(ROIdata,'PCA')
    load([movementdir,files(filei).name]);
    for videoi=1:length(rawvideo)
        %         if videoi==1
        %             bigvideo=rawvideo(videoi).vid;
        %         else
        %             bigvideo=cat(3,bigvideo,rawvideo(videoi).vid);
        %         end
        bigvideo=rawvideo(videoi).vid;
        %%
        for ROIi=1:length(ROIdata)
            if isempty(ROIdata(ROIi).mask)
                disp(['no ROI selected for ROI:',ROIdata(ROIi).ROIname,' in file ',files(filei).name]);
            elseif overwrite==1 | ~isfield(ROIdata,'PCA') | ~isfield(ROIdata(ROIi).PCA,'mask') | length(ROIdata(ROIi).PCA)<videoi | isempty(ROIdata(ROIi).PCA(videoi).mask) | ROIdata(ROIi).PCA(videoi).mask~=ROIdata(ROIi).mask | ~isfield(ROIdata(ROIi).PCA(videoi),'uniquePC')
                ROIdata(ROIi).PCA(videoi).mask=ROIdata(ROIi).mask;
                bigmask=repmat(ROIdata(ROIi).mask,1,1,size(bigvideo,3));
                ROIvideo=bigvideo;
                ROIvideo(~bigmask)=0;
                %                 tic
                linearROIvideo=reshape(ROIvideo,size(ROIvideo,1)*size(ROIvideo,2),size(ROIvideo,3));
                %%
                borders=[];
                whattoprincomp=((double(linearROIvideo')));
                whattoprincomp=diff(zscore(whattoprincomp));
                %%
                % %                 borders(1,:) = prctile(whattoprincomp,2);
                % %                 borders(2,:) = prctile(whattoprincomp,98);
                % %                 borders(3,:) = diff(borders);
                % % %                 borders=borders;
                % %                 whattoprincomp=whattoprincomp-repmat(borders(1,:),size(whattoprincomp,1),1);
                % %                 whattoprincomp=whattoprincomp./repmat(borders(3,:),size(whattoprincomp,1),1);
                % %                 whattoprincomp(whattoprincomp>1)=1;
                % %                 whattoprincomp(whattoprincomp<0)=0;
                % %                 whattoprincomp=diff(whattoprincomp);
                % %                 whattoprincompold=whattoprincomp;
                % %                 rowstodel=find(isnan(whattoprincomp(end,:)));
                
                rowstodel=find(sum(whattoprincomp)==0);
                whattoprincomp(:,rowstodel)=[];
                
                %%
                [coeff,score,latent]=princomp(whattoprincomp);
                ROIdata(ROIi).PCA(videoi).PC=zeros(size(score,1),10);
                eddig=min(10,size(score,2));
                ROIdata(ROIi).PCA(videoi).PC(:,1:eddig)=score(:,1:eddig);
                ROIdata(ROIi).PCA(videoi).latent=latent;
                ROIdata(ROIi).PCA(videoi).cumweigths=cumsum(latent)/sum(latent);
                ROIdata(ROIi).PCA(videoi).weigths=(latent)/sum(latent);
                %                 toc
                %         %%
                %         scorelayout=[];
                %         close all
                %         princompstoshow=600;
                %         figure(10)
                %         clf
                %         for i=1:princompstoshow
                %             %                 subplot(princompstoshow,2,(i-1)*2+1)
                %             subplot(1,2,1)
                %             plot((score(:,i)))
                %             %                 subplot(princompstoshow,2,(i-1)*2+2)
                %             subplot(1,2,2)
                % %             scorelayout(:,:,i)=reshape(coeff(:,i),size(bigvideo,1),size(bigvideo,2));
                % %             imagesc(scorelayout(:,:,i))
                %             %                 caxis([-.04 .04])
                %             pause
                %         end
                %%
                if strcmp(ROIdata(ROIi).ROIname,'Body')
                    powervals=[];
                    for pci=1:10
                        parameters=struct;
                        parameters.min=1; %minimal frequency for decomposition
                        parameters.max=5;% - maximal frequency for decomposition
                        parameters.step=.1; %- frequency stepsize
                        parameters.scale=1;% - scale type, can be linear (1) or logarithmic (2)
                        parameters.wavenumber=7;% - number of waves in wavelet
                        parameters.waveletlength=10;
                        parameters.addtaper=1;
                        parameters.taperlength=5;
                        parameters.interval=mode(diff(videodata(videoi).time));
                        [powerMatrix, frequencyVector]  = calculateTFDPower(ROIdata(ROIi).PCA(videoi).PC(:,pci), parameters);
                        powervals(pci)=median(max(powerMatrix));
%                         figure(1)
%                         clf
%                         subplot(3,1,1)
%                         imagesc(1:length(ROIdata(ROIi).PCA(videoi).PC(:,pci)),frequencyVector,powerMatrix);
%                         set(gca,'YDir','normal');
%                         subplot(3,1,2)
%                         plot(max(powerMatrix));
%                         hold on
%                         plot(max(powerMatrix)*0+median(max(powerMatrix)),'r-')
%                         subplot(3,1,3)
%                         plot(ROIdata(ROIi).PCA(videoi).PC(:,pci))
%                         pause
                    end
                    [~,breathindex]=max(powervals);
                    ROIdata(ROIi).PCA(videoi).breath=ROIdata(ROIi).PCA(videoi).PC(:,breathindex);
                end
            end
        end
        %% PCA based PC selection
        ennyiPCkell=3;
        PCs=struct;
        for ROIi=1:length(ROIdata)
            PCs(ROIi).pc=ROIdata(ROIi).PCA(videoi).PC(:,1:ennyiPCkell);
        end
        whattoprincomp=zscore([PCs.pc]);
        [coeff,score,latent]=princomp(whattoprincomp);
        osszidx=1:length(ROIdata)*ennyiPCkell;
        for ROIi=4:length(ROIdata)
            idxmost=(ROIi-1)*ennyiPCkell+1:ennyiPCkell+(ROIi-1)*ennyiPCkell;
            idxrest=osszidx;
            idxrest(idxmost)=[];
            idxbckground=[-ennyiPCkell+1:0]+ennyiPCkell*length(ROIdata);
            [bestcompcoeff,bestcompidx]=max(abs(coeff(idxmost,:)),[],2);
            
            [correspondingbackgroundcoeff,corrbcgidx]=max(abs(coeff(idxbckground,bestcompidx)));
            [correspondingrestcoeff,~]=max(abs(coeff(idxrest,bestcompidx)));
            
            coeffbackgroundratio=bestcompcoeff./correspondingbackgroundcoeff';
            coeffrestratio=bestcompcoeff./correspondingrestcoeff';
            
            [~,idx]=max(coeffbackgroundratio);
            ROIdata(ROIi).PCA(videoi).bestPC=ROIdata(ROIi).PCA(videoi).PC(:,idx);
            [~,idx]=max(coeffrestratio);
            ROIdata(ROIi).PCA(videoi).uniquePC=ROIdata(ROIi).PCA(videoi).PC(:,idx);
            
%             figure(1)
%             clf
%             for princompi=1:ennyiPCkell
%                 subplot(ennyiPCkell,1,princompi)
%                 plot(whattoprincomp(:,idxmost(princompi)))
%                 hold on
%                 plot(whattoprincomp(:,idxbckground(corrbcgidx(princompi))),'r-')
%             end
            
            %%% itt kiválasztom azt a principális komponenst, ami egy adott
            %%% ROI-ra jellemző, de a háttérre nem
            %             return
        end
        
        %%
        %         for i=1:100
        %
        %                 subplot(2,2,1)
        %                 plot((score(:,i)))
        %                xlabel('Frames')
        %                 ylabel('master PCA score')
        %
        %                 subplot(2,2,2)
        %                 plot([1:ennyiPCkell*length(ROIdata)]/ennyiPCkell+1,(coeff(:,i)))
        %                 xlabel('ROIs')
        %                 ylabel('PCA coefficient')
        %
        %                 [~,idx]=max(abs(coeff(:,i)));
        %                 subplot(2,2,3)
        %                 plot(whattoprincomp(:,idx));
        %
        %                 xlabel('Frames')
        %                 ylabel('original PCA score')
        %                 subplot(2,2,4)
        %                 imagesc(abs(coeff'))
        %                 set(gca,'YDir','normal')
        %                 pause
        %             end
    end
    disp(['PCA for ', files(filei).name, 'is done'])
    save([ROIdir,files(filei).name],'ROIdata');
    else
        disp(['no ROIS selected for ', files(filei).name])
    end 
end