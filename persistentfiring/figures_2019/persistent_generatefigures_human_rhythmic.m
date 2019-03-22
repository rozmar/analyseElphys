%% Figure 1
%% HUMAN nonrhythmic
ID='1603021rm_5_3,4,5_4';
expname=[ID,'HUMAN_nonrhythmic_persistent'];
xlsidx=find(strcmp({xlsdata.ID},ID));
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
persistent_generatefigures_plotselectedtimeinterval(xlsidx,dirs,xlsdata,valtozok,expname,additionaldata);

%% HUMAN rhythmic 
ID='1702021rm_6_1_3';
xlsidx=find(strcmp({xlsdata.ID},ID));
expname=[ID,'_HUMAN_rhythmic_persistent'];
valtozok_sampletrace=struct;
valtozok_sampletrace.debugmode=0;
valtozok_sampletrace.zerotime=65124.6;
valtozok_sampletrace.xlimits=[-5 38.5]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.xlimitstoexclude=[0,0]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.xlimitsblowup=[-4,-.5;20,23.5]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.ylimitsblowup=zeros(size(valtozok_sampletrace.xlimitsblowup));
valtozok_sampletrace.ylimits=[-73,31];%
% valtozok_sampletrace.ylimitscurr=[0 0];%[-150,250];
valtozok_sampletrace.isiYlimits=[0 2];
valtozok_sampletrace.freqYlimits=[.1 500];
valtozok_sampletrace.freqYscale=[0.1, 1, 10, 100];
valtozok_sampletrace.cutofffreq=8500;
valtozok_sampletrace.highlightaxonalspikes=1;
valtozok_sampletrace.highlightaxonalspikes_timeback=.001;
valtozok_sampletrace.highlightaxonalspikes_timeforward=.010;
valtozok_sampletrace.xcm=17;
valtozok_sampletrace.xcm_blowup=8;
valtozok_sampletrace.ycm=2;
valtozok_sampletrace.ycm_current=1;
valtozok_sampletrace.voltagelinewidth=.5;
valtozok_sampletrace.currentlinewidth=1;
valtozok_sampletrace.fontsize=8;
valtozok_sampletrace.fonttype='Helvetica';
valtozok_sampletrace.axeswidth=1;
persistent_generatefigures_plotselectedtimeinterval(xlsidx,dirs,xlsdata,valtozok_sampletrace,expname);

%% APwaves generation
ID='1702021rm_6_1_3';
xlsidx=find(strcmp({xlsdata.ID},ID));
valtozok_APwaves=struct;
valtozok_APwaves.debugmode=0;
valtozok_APwaves.zerotime=65124.6;
valtozok_APwaves.xlimits=valtozok_sampletrace.xlimits;%[-5 5]+valtozok_APwaves.zerotime;
valtozok_APwaves.filter='gauss';% 'gauss' - gaussian filter, 'boxcar' - moving average
valtozok_APwaves.movingt=20;%microseconds for filtering: sigma in case of gaussian filter, width in case of boxcar
valtozok_APwaves.zerosettime='threshold'; %threshold, apmaximum
valtozok_APwaves.timeback=.0005; %ms back from zero time
valtozok_APwaves.timeforward=.002; %ms forward from zero time
APwaves=persistent_generatefigures_generate_APwaves(dirs,xlsdata,xlsidx,valtozok_APwaves);
needed=(valtozok_sampletrace.xlimitsblowup(1,1)<=[APwaves.maxtime] & valtozok_sampletrace.xlimitsblowup(2,1)>=[APwaves.maxtime]) | (valtozok_sampletrace.xlimitsblowup(1,2)<=[APwaves.maxtime] & valtozok_sampletrace.xlimitsblowup(2,2)>=[APwaves.maxtime]);
APwaves=APwaves(needed);
%%
valtozok_APwaves_plot=struct;
valtozok_APwaves_plot.markersize=5;
valtozok_APwaves_plot.linewidth=.5;
valtozok_APwaves_plot.highlight.Window=[0, 0] + valtozok_APwaves.zerotime;
% valtozok_APwaves_plot.highlight.idx=[32,100];% indexes to highlight
aptimestohighlight=[-3;21.5];
[~,valtozok_APwaves_plot.highlight.idx]=min(abs(repmat([APwaves.maxtime]-valtozok_APwaves.zerotime,length(aptimestohighlight),1)-repmat(aptimestohighlight,1,length([APwaves.maxtime]))),[],2);
valtozok_APwaves_plot.highlight.linewidth=1;
valtozok_APwaves_plot.xcm=4;
valtozok_APwaves_plot.ycm=4;
valtozok_APwaves_plot.fontsize=8;
valtozok_APwaves_plot.fonttype='Helvetica';
valtozok_APwaves_plot.axeswidth=1;
persistent_generatefigures_plot_APwaves(xlsidx,dirs,xlsdata,valtozok_APwaves_plot,expname,APwaves)%(APwaves,valtozok_APwaves_plot)