%% HUMAN nonrhythmic

ID='1603021rm_5_3,4,5_4';
expname=[ID,'HUMAN_nonrhythmic_persistent'];
xlsidx=find(strcmp({xlsdata.ID},ID));
load([dirs.eventdir,'sorted/',ID,'.mat'],'eventdata');

valtozok=struct;
valtozok.debugmode=0;
valtozok.zerotime=61484.85;%61391.4+92+1.45;
valtozok.xlimits=[-5 38.5]+valtozok.zerotime;
valtozok.xlimitstoexclude=[0,0]+valtozok.zerotime;
valtozok.xlimitsblowup=[0,0]+valtozok.zerotime;
valtozok.ylimitsblowup=[0,0];
valtozok.ylimits=[0,0];%
% valtozok.ylimitscurr=[0 0];%[-150,250];
valtozok.isiYlimits=[0 2];
valtozok.freqYlimits=[.1 500];
valtozok.freqYscale=[0.1, 1, 10, 100];
valtozok.cutofffreq=8500;
valtozok.highlightaxonalspikes=1;
valtozok.highlightaxonalspikes_timeback=.001;
valtozok.highlightaxonalspikes_timeforward=.010;
valtozok.xcm=17;
valtozok.ycm=2;
valtozok.ycm_current=1;
valtozok.fontsize=8;
valtozok.fonttype='Helvetica';
additionaldata=struct;
additionaldata.eventdata=eventdata;
persistent_plotselectedtimeinterval(xlsidx,dirs,xlsdata,valtozok,expname,additionaldata);


%% HUMAN rhythmic
ID='1702021rm_6_1_3';
xlsidx=find(strcmp({xlsdata.ID},ID));
expname=[ID,'_HUMAN_rhythmic_persistent'];

load([dirs.eventdir,'sorted/',ID,'.mat'],'eventdata');

valtozok=struct;
valtozok.debugmode=0;
valtozok.zerotime=65124.6;
valtozok.xlimits=[-5 38.5]+valtozok.zerotime;
valtozok.xlimitstoexclude=[0,0]+valtozok.zerotime;
valtozok.xlimitsblowup=[-4,-.5;20,23.5]+valtozok.zerotime;
valtozok.ylimitsblowup=zeros(size(valtozok.xlimitsblowup));
valtozok.ylimits=[-73,31];%
% valtozok.ylimitscurr=[0 0];%[-150,250];
valtozok.isiYlimits=[0 2];
valtozok.freqYlimits=[.1 500];
valtozok.freqYscale=[0.1, 1, 10, 100];
valtozok.cutofffreq=8500;
valtozok.highlightaxonalspikes=1;
valtozok.highlightaxonalspikes_timeback=.001;
valtozok.highlightaxonalspikes_timeforward=.010;
valtozok.xcm=17;
valtozok.xcm_blowup=8;
valtozok.ycm=2;
valtozok.ycm_current=1;
valtozok.fontsize=8;
valtozok.fonttype='Helvetica';
valtozok.axeswidth=1;
additionaldata=struct;
additionaldata.eventdata=eventdata;
persistent_plotselectedtimeinterval(xlsidx,dirs,xlsdata,valtozok,expname,additionaldata);

%% APwaves generation
ID='1702021rm_6_1_3';
xlsidx=find(strcmp({xlsdata.ID},ID));
valtozok=struct;
valtozok.debugmode=0;
valtozok.zerotime=65124.6;
valtozok.xlimits=[-5 5]+valtozok.zerotime;
valtozok.filter='gauss';% 'gauss' - gaussian filter, 'boxcar' - moving average
valtozok.movingt=15;%microseconds for filtering: sigma in case of gaussian filter, width in case of boxcar
valtozok.zerosettime='threshold'; %threshold, apmaximum
valtozok.timeback=.0005; %ms back from zero time
valtozok.timeforward=.0015; %ms forward from zero time
APwaves=persistent_generatefigures_generate_APwaves(dirs,xlsdata,xlsidx,valtozok);

%%
valtozok_plot=struct;
valtozok_plot.markersize=10;
valtozok_plot.linewidth=1;
valtozok_plot.emphasize.Window=[0, 0] + valtozok.zerotime;
valtozok_plot.emphasize.idx=[32,100];% indexes to highlight
valtozok_plot.emphasize.linewidth=3;

valtozok_plot.axeswidth=1;
% valtozok_plot.timeborders=valtozok.timeborders;
persistent_generatefigures_plot_APwaves(APwaves,valtozok_plot)