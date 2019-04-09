function persistent_generatefigures_plot_PSD_for_timeinterval(xlsidx,dirs,xlsdata,valtozok,expname)
parameters=valtozok.PSDparameters;
if isfield(valtozok,'xcm')
    xcm=valtozok.xcm;
else
    xcm=20;
end
if isfield(valtozok,'xcm_blowup')
    xcm_blowup=valtozok.xcm_blowup;
else
    xcm_blowup=xcm;
end
if isfield(valtozok,'ycm')
    ycm=valtozok.ycm;
else
    ycm=10;
end

if isfield(valtozok,'dpi')
    dpi=valtozok.dpi;
else
    dpi=900;
end
if isfield(valtozok,'fontsize')
    betumeret=valtozok.fontsize;
else
    betumeret=14;
end
if isfield(valtozok,'fonttype')
    betutipus=valtozok.fonttype;
else
    betutipus='Helvetica';
end
if isfield(valtozok,'axeswidth')
    axesvastagsag=valtozok.axeswidth;
else
    axesvastagsag=1;
end
if isfield(valtozok,'caxvals')
    caxvals=valtozok.caxvals;
else
    caxvals=[];
end
if ~isfield(valtozok,'freqticks')
    valtozok.freqticks=[];
end

renderer='painters';


fname=[xlsdata(xlsidx).ID,'.mat'];

load([dirs.bridgeddir,fname],'bridgeddata','stimdata');

neededidx=([bridgeddata.realtime]>valtozok.xlimits(1) &[bridgeddata.realtime]<valtozok.xlimits(2));
if isfield(valtozok,'xlimitstoexclude')
    for i=1:size(valtozok.xlimitstoexclude,1)
        neededidx=neededidx &~([bridgeddata.realtime]>valtozok.xlimitstoexclude(i,1) &[bridgeddata.realtime]<valtozok.xlimitstoexclude(i,2));
    end
end
neededidx=find(neededidx);
if neededidx(1)~=1
    neededidx=[neededidx(1)-1,neededidx];
end
%%
bridgedsweeps=struct;
for sweep=1:length(neededidx)
    
    bridgedsweeps(sweep).starttime=bridgeddata(neededidx(sweep)).realtime;
    
    
    
    bridgedsweeps(sweep).sweepnum=neededidx(sweep);
    
    
    si=bridgeddata(neededidx(sweep)).si;
    cutoffreq=parameters.max*10;
    [b,a]=butter(1,cutoffreq/(1/si)/2,'low');
    y = filtfilt(b,a,bridgeddata(neededidx(sweep)).y);
    y=downsample(y,round((1/cutoffreq)/si));
    bridgedsweeps(sweep).y=y;
    newsi=1/cutoffreq;
    bridgedsweeps(sweep).si=newsi;   
    bridgedsweeps(sweep).x=[bridgeddata(neededidx(sweep)).realtime:newsi:bridgeddata(neededidx(sweep)).realtime+(length(bridgedsweeps(sweep).y)-1)*newsi]-valtozok.zerotime;
    
     
end
%% 
clc
si=1/cutoffreq;
y_cont=bridgedsweeps(1).y;
x_cont=bridgedsweeps(1).x;

if length(bridgedsweeps)>1
    for sweepi=2:length(bridgedsweeps)
        timedifi=bridgedsweeps(sweepi).x(1)-x_cont(end)+1*si;
        stepdifi=round(timedifi/si);
        x_cont=[x_cont,x_cont(end)+[1:stepdifi]*si,bridgedsweeps(sweepi).x];
        y_cont=[y_cont,y_cont(end)+[1:stepdifi]/stepdifi*(bridgedsweeps(sweepi).y(1)-y_cont(end)),bridgedsweeps(sweepi).y];
        
        a=[1:stepdifi]/stepdifi/timedifi;
        mode(diff(a)) 
    end
end
parameters.interval=newsi;
[powerMatrix, frequencyVector]  = calculateTFDPower(y_cont', parameters);
figure(3)
clf
imagesc(x_cont,frequencyVector,powerMatrix)
box off
set(gca,'YDir','normal');
colormap linspecer
ylabel('Frequency (Hz)')
xlabel('Time (s)')
xlim(valtozok.xlimits-valtozok.zerotime)
if ~isempty(caxvals)
    caxis(caxvals);
end





set(gca,'LineWidth',axesvastagsag,'FontSize',betumeret,'Fontname',betutipus,'Units','normalized','Position',[.25 .25 .5 .5])
set(findobj(gcf,'type','text'),'fontsize',betumeret,'Fontname',betutipus)
set(gcf,'PaperUnits','centimeters','PaperPositionMode','manual','PaperSize',[xcm/.5 ycm/.5]+2,'PaperPosition',[2 2 xcm/.5 ycm/.5])



if parameters.scale==2 & ~isempty(valtozok.freqticks)
    ylimnow=get(gca,'Ylim');
    linticks = linspace(ylimnow(1),ylimnow(end),1000);
    logticks = logspace(log10(ylimnow(1)),log10(ylimnow(end)),1000);
    
    for i=1:length(valtozok.freqticks)
        [~,tickidx(i)]=min(abs(logticks-valtozok.freqticks(i)));
    end
    set(gca,'Ytick',linticks(tickidx))
    set(gca,'Yticklabel',valtozok.freqticks)
elseif parameters.scale==2
    set(gca, 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
    ytickidxs=[];
    yticknow=get(gca,'YTick');
    ylimnow=get(gca,'Ylim');
    
    linearyticks = linspace(ylimnow(1),ylimnow(end),length(frequencyVector));
    % linearyticks = frequencyVector;
    for ticki=1:length(yticknow)
        [~,ytickidxs(ticki)]=min(abs(linearyticks-yticknow(ticki)));
    end
    yticklabelnow=frequencyVector(ytickidxs);
    set(gca,'Yticklabel',round(yticklabelnow*100)/100)
end

%%




%%


%%
set(gcf, 'Renderer', renderer);
% saveas(gcf,[dirs.figuresdir,expname,'_voltage.pdf'])
print(gcf,[dirs.figuresdir,expname,'_PSD.jpg'],'-djpeg',['-r',num2str(dpi)])