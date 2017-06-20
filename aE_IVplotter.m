%%
close all
locations=marcicucca_locations;
% clear all
xbar=.1; %s
ybar=.02; %V
setup='2P3DAO';
file='1310302rm';
gsc='g1_s1_c3';
sweeps=[1,15];
overwrite=1;
for xlsi=1:length(xlsdata)
    figure(1)
    clf
    
    if xlsdata(xlsi).field==0 &&  xlsdata(xlsi).juxta==0
        setup=xlsdata(xlsi).setup;
        file=xlsdata(xlsi).HEKAfname;
        IVstring=xlsdata(xlsi).G_S_C;
        
        gs=[];
        ss=[];
        cs=[];
        hyps=strfind(IVstring,'_');
        gs=str2num(IVstring(1:hyps(1)-1));
        ss=str2num(IVstring(hyps(1)+1:hyps(2)-1));
        cs=str2num(IVstring(hyps(2)+1:end));
        if ischar(xlsdata(xlsi).Ivsweeps)
            secondsweeps=str2num(xlsdata(xlsi).Ivsweeps);
        else
            secondsweeps=xlsdata(xlsi).Ivsweeps;
        end
        for si=1:length(ss)
            gsc=['g',num2str(gs),'_s',num2str(ss(si)),'_c',num2str(cs)];
            a=dir([dirs.figuresdir,'IV_',file,'_',gsc,'.pdf']);
            if (isempty(a) | overwrite==1)
                sweeps=[1,secondsweeps(si)];
                
                filter=[2000,15000];
                
                aa=dir([locations.tgtardir,'MATLABdata/IV/',setup,'/',file,'.mat']);
                if ~isempty(aa)
                    load([locations.tgtardir,'MATLABdata/IV/',setup,'/',file]);
                    iv=iv.(gsc);
                    figure(1)
                    clf
                    hold on;
                    for i=1:length(sweeps)
                        si=mode(diff(iv.time));
                        
                        v=iv.(['v',num2str(sweeps(i))]);
                        if  filter(min(length(filter),i))>0
                            [b,a]=butter(1,filter(i)/(1/mode(diff(iv.time)))/2,'low');
                            v=filtfilt(b,a,v);
                        end
                        plot(iv.time,v,'k-','LineWidth',2);
                        
                    end
                    axis tight
                    xlimits=get(gca,'Xlim');
                    ylimits=get(gca,'Ylim');
                    yval=mean(ylimits);
                    xval=max(xlimits);
                    plot([xval,xval+xbar],[yval,yval],'k-','LineWIdth',4)
                    plot([xval,xval],[yval,yval+ybar],'k-','LineWIdth',4)
                    text(xval+ybar,yval+ybar,['\bf',num2str(ybar*1000),' mV'])
                    text(xval+ybar,yval+ybar/2,['\bf',num2str(xbar*1000),' ms'])
                    %     legend('40 mV')
                    axis off
                    axis tight
                    %                 plot2svg([dirs.figuresdir,'IV_',file,'_',gsc,'.svg'],gcf);
                    print([dirs.figuresdir,'IV_',file,'_',gsc,'.pdf'],'-dpdf')
                end
            end
        end
    end
    % 
end