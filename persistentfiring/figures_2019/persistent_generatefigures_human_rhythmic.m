%% Figure 1
% These scripts export figures in .jpg format which can be inserted and
% further organized in Inkscape. The figures are saved in the
% dirs.figuresdir folder.
% see also persistent_generatefigures_plotselectedtimeinterval,
% persistent_generatefigures_generate_APwaves, persistent_generatefigures_plot_APwaves,
% persistent_generatefigures_plot_PSD_for_timeinterval
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
valtozok.markersize=1;
valtozok.fontsize=8;
valtozok.fonttype='Helvetica';
valtozok.axeswidth=.5;
valtozok.axis.voltage_y=true;
valtozok.axis.voltage_x=false;
valtozok.axis.current_y=true;
valtozok.axis.current_x=false;
valtozok.axis.freq_y=true;
valtozok.axis.freq_x=true;
persistent_generatefigures_plotselectedtimeinterval(xlsidx,dirs,xlsdata,valtozok,expname);

%% HUMAN sporadic
ID='1602261rm_5_1_4';
expname=[ID,'HUMAN_sporadic_persistent'];
xlsidx=find(strcmp({xlsdata.ID},ID));
valtozok=struct;
valtozok.debugmode=0;
valtozok.zerotime=73271;%61391.4+92+1.45;
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
valtozok.highlightaxonalspikes_timeforward=.015;
valtozok.xcm=17;
valtozok.ycm=2;
valtozok.ycm_current=1;
valtozok.markersize=1;
valtozok.fontsize=8;
valtozok.fonttype='Helvetica';
valtozok.axeswidth=.5;
valtozok.axis.voltage_y=true;
valtozok.axis.voltage_x=false;
valtozok.axis.current_y=true;
valtozok.axis.current_x=false;
valtozok.axis.freq_y=true;
valtozok.axis.freq_x=true;
persistent_generatefigures_plotselectedtimeinterval(xlsidx,dirs,xlsdata,valtozok,expname);

%% HUMAN rhythmic 
ID='1702021rm_6_1_3';
xlsidx=find(strcmp({xlsdata.ID},ID));
expname=[ID,'_HUMAN_rhythmic_persistent'];
valtozok_sampletrace=struct;
valtozok_sampletrace.debugmode=0;
valtozok_sampletrace.zerotime=65124.6;
valtozok_sampletrace.xlimits=[-5 38.5]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.xlimitstoexclude=[0,0]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.xlimitsblowup=[-4,-1.5;21,23.5]+valtozok_sampletrace.zerotime;%[-4,-.5;20,23.5]+valtozok_sampletrace.zerotime;%
valtozok_sampletrace.ylimitsblowup=zeros(size(valtozok_sampletrace.xlimitsblowup));
valtozok_sampletrace.ylimits=[-73,35];%
% valtozok_sampletrace.ylimitscurr=[0 0];%[-150,250];
valtozok_sampletrace.isiYlimits=[0 2];
valtozok_sampletrace.freqYlimits=[.1 500];
valtozok_sampletrace.freqYscale=[0.1, 1, 10, 100];
valtozok_sampletrace.cutofffreq=8500;
valtozok_sampletrace.highlightaxonalspikes=1;
valtozok_sampletrace.highlightaxonalspikes_timeback=.001;
valtozok_sampletrace.highlightaxonalspikes_timeforward=.020;
valtozok_sampletrace.xcm=17;
valtozok_sampletrace.xcm_blowup=8;
valtozok_sampletrace.ycm=2;
valtozok_sampletrace.ycm_current=.5;
valtozok_sampletrace.voltagelinewidth=.5;
valtozok_sampletrace.currentlinewidth=.5;
valtozok_sampletrace.markersize=1;
valtozok_sampletrace.fontsize=8;
valtozok_sampletrace.fonttype='Helvetica';
valtozok_sampletrace.axeswidth=.5;
valtozok_sampletrace.axis.voltage_y=true;
valtozok_sampletrace.axis.voltage_x=false;
valtozok_sampletrace.axis.current_y=false;
valtozok_sampletrace.axis.current_x=false;
valtozok_sampletrace.axis.freq_y=true;
valtozok_sampletrace.axis.freq_x=true;

persistent_generatefigures_plotselectedtimeinterval(xlsidx,dirs,xlsdata,valtozok_sampletrace,expname);

%% APwaves generation - first run sampletrace extraction for time boundaries
ID='1702021rm_6_1_3';
xlsidx=find(strcmp({xlsdata.ID},ID));
valtozok_APwaves=struct;
valtozok_APwaves.debugmode=0;
valtozok_APwaves.zerotime=65124.6;
valtozok_APwaves.xlimits=valtozok_sampletrace.xlimits;%[-5 5]+valtozok_APwaves.zerotime;
valtozok_APwaves.filter='gauss';% 'gauss' - gaussian filter, 'boxcar' - moving average
valtozok_APwaves.movingt=20;%microseconds for filtering: sigma in case of gaussian filter, width in case of boxcar
valtozok_APwaves.zerosettime='threshold'; %threshold, apmaximum
valtozok_APwaves.timeback=.001; %ms back from zero time
valtozok_APwaves.timeforward=.002; %ms forward from zero time
APwaves=persistent_generatefigures_generate_APwaves(dirs,xlsdata,xlsidx,valtozok_APwaves);
needed=(valtozok_sampletrace.xlimitsblowup(1,1)<=[APwaves.maxtime] & valtozok_sampletrace.xlimitsblowup(1,2)>=[APwaves.maxtime]) | (valtozok_sampletrace.xlimitsblowup(2,1)<=[APwaves.maxtime] & valtozok_sampletrace.xlimitsblowup(2,2)>=[APwaves.maxtime]);
APwaves=APwaves(needed);
% APwaves([APwaves.stimulated]&[APwaves.axonalAP])=[]; %% axonal aps from stimulated ss
%
valtozok_APwaves_plot=struct;
valtozok_APwaves_plot.markersize=5;
valtozok_APwaves_plot.linewidth=.5;
valtozok_APwaves_plot.highlight.Window=[0, 0] + valtozok_APwaves.zerotime;
% valtozok_APwaves_plot.highlight.idx=[32,100];% indexes to highlight
aptimestohighlight=[-3;21.5];
[~,valtozok_APwaves_plot.highlight.idx]=min(abs(repmat([APwaves.maxtime]-valtozok_APwaves.zerotime,length(aptimestohighlight),1)-repmat(aptimestohighlight,1,length([APwaves.maxtime]))),[],2);
valtozok_APwaves_plot.highlight.linewidth=1;
valtozok_APwaves_plot.xcm=3;
valtozok_APwaves_plot.ycm=3;
valtozok_APwaves_plot.xcm_scatterplot=5;
valtozok_APwaves_plot.ycm_scatterplot=5;
valtozok_APwaves_plot.fontsize=8;
valtozok_APwaves_plot.fonttype='Helvetica';
valtozok_APwaves_plot.axeswidth=.5;
valtozok_APwaves_plot.separateapwaves=1;
% valtozok_sampletrace.axis.voltage_y=true;
% valtozok_sampletrace.axis.voltage_x=false;
% valtozok_sampletrace.axis.current_y=true;
% valtozok_sampletrace.axis.current_x=false;
% valtozok_sampletrace.axis.freq_y=true;
% valtozok_sampletrace.axis.freq_x=true;
persistent_generatefigures_plot_APwaves(xlsidx,dirs,xlsdata,valtozok_APwaves_plot,expname,APwaves)%(APwaves,valtozok_APwaves_plot)

%% Figure 2
%% pharmacology - NBQX, GBZ, APV, CGP - sample trace
ID='1411283rm_2_1,4_4';
xlsidx=find(strcmp({xlsdata.ID},ID));
expname=[ID,'_HUMAN_rhythmic_persistent_NBQX_GBZ_APV_CGP'];
valtozok_sampletrace=struct;
valtozok_sampletrace.debugmode=0;
valtozok_sampletrace.zerotime=12650;
valtozok_sampletrace.xlimits=[0 1350]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.ylimits=[-72 31];%
valtozok_sampletrace.xlimitstoexclude=[0,0]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.xlimitsblowup=[38.7, 48.7;678.7,688.7;1214, 1224]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.ylimitsblowup=zeros(size(valtozok_sampletrace.xlimitsblowup));
valtozok_sampletrace.ylimitsblowup=[-72,24;-72,24;-72,24];

% valtozok_sampletrace.ylimitscurr=[0 0];%[-150,250];
valtozok_sampletrace.isiYlimits=[0 2];
valtozok_sampletrace.freqYlimits=[.1 500];
valtozok_sampletrace.freqYscale=[0.1, 1, 10, 100];
valtozok_sampletrace.cutofffreq=8500;
valtozok_sampletrace.highlightaxonalspikes=1;
valtozok_sampletrace.highlightaxonalspikes_timeback=.001;
valtozok_sampletrace.highlightaxonalspikes_timeforward=.020;
valtozok_sampletrace.drugwashin.drugwashinline_ystart=26;
valtozok_sampletrace.drugwashin.drugwashinline_ystep=5;
valtozok_sampletrace.drugwashin.drugwashinline_linewidth=2;

valtozok_sampletrace.xcm=17;
valtozok_sampletrace.xcm_blowup=5;
valtozok_sampletrace.ycm=2;
valtozok_sampletrace.ycm_current=1;
valtozok_sampletrace.voltagelinewidth=.25;
valtozok_sampletrace.currentlinewidth=1;
valtozok_sampletrace.markersize=1;
valtozok_sampletrace.fontsize=8;
valtozok_sampletrace.fonttype='Helvetica';
valtozok_sampletrace.axeswidth=.5;
valtozok_sampletrace.axis.voltage_y=true;
valtozok_sampletrace.axis.voltage_x=false;
valtozok_sampletrace.axis.current_y=true;
valtozok_sampletrace.axis.current_x=false;
valtozok_sampletrace.axis.freq_y=true;
valtozok_sampletrace.axis.freq_x=true;

persistent_generatefigures_plotselectedtimeinterval(xlsidx,dirs,xlsdata,valtozok_sampletrace,expname);

%% pharmacology - NBQX, GBZ, APV, CGP - PSD plot
ID='1411283rm_2_1,4_4';
xlsidx=find(strcmp({xlsdata.ID},ID));
expname=[ID,'_HUMAN_rhythmic_persistent_NBQX_GBZ_APV_CGP'];
valtozok_PSD=struct;
valtozok_PSD.debugmode=0;
valtozok_PSD.zerotime=12650;
valtozok_PSD.xlimits=[0 1350]+valtozok_PSD.zerotime;
valtozok_PSD.xlimitstoexclude=[0,0]+valtozok_PSD.zerotime;

valtozok_PSD.xcm=17;
valtozok_PSD.xcm_blowup=5;
valtozok_PSD.ycm=2;
valtozok_PSD.fontsize=8;
valtozok_PSD.fonttype='Helvetica';
valtozok_PSD.axeswidth=.5;
valtozok_PSD.caxvals=[0 .1];
valtozok_PSD.freqticks=[0.3,1,3,9];
% valtozok_PSD.freqticks=[0.5,1,2,4];

valtozok_PSD.PSDparameters=struct;
valtozok_PSD.PSDparameters.min=.25; %minimal frequency for decomposition
valtozok_PSD.PSDparameters.max=10;% - maximal frequency for decomposition
valtozok_PSD.PSDparameters.step=.01; %- frequency stepsize
valtozok_PSD.PSDparameters.scale=2;% - scale type, can be linear (1) or logarithmic (2)
valtozok_PSD.PSDparameters.wavenumber=9;% - number of waves in wavelet
valtozok_PSD.PSDparameters.waveletlength=30;
valtozok_PSD.PSDparameters.addtaper=1;
valtozok_PSD.PSDparameters.taperlength=30;

persistent_generatefigures_plot_PSD_for_timeinterval(xlsidx,dirs,xlsdata,valtozok_PSD,expname);

%% voltage dependence - sample trace

ID='1411283rm_2_1,4_4';
xlsidx=find(strcmp({xlsdata.ID},ID));
expname=[ID,'_HUMAN_rhythmic_persistent_voltage_dependence'];
valtozok_sampletrace=struct;
valtozok_sampletrace.debugmode=0;
valtozok_sampletrace.zerotime=11838.5;
valtozok_sampletrace.xlimits=[-5 90]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.ylimits=[-101 28];%
valtozok_sampletrace.xlimitstoexclude=[0,0]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.xlimitsblowup=[25, 30;51,56;65, 70;87,92]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.ylimitsblowup=zeros(size(valtozok_sampletrace.xlimitsblowup));
% valtozok_sampletrace.ylimitsblowup=[-72,24;-72,24;-72,24];

% valtozok_sampletrace.ylimitscurr=[0 0];%[-150,250];
valtozok_sampletrace.isiYlimits=[0 2];
valtozok_sampletrace.freqYlimits=[.1 500];
valtozok_sampletrace.freqYscale=[1, 10, 100];
valtozok_sampletrace.cutofffreq=8500;
valtozok_sampletrace.highlightaxonalspikes=1;
valtozok_sampletrace.highlightaxonalspikes_timeback=.001;
valtozok_sampletrace.highlightaxonalspikes_timeforward=.020;
valtozok_sampletrace.drugwashin.drugwashinline_ystart=26;
valtozok_sampletrace.drugwashin.drugwashinline_ystep=5;
valtozok_sampletrace.drugwashin.drugwashinline_linewidth=2;

valtozok_sampletrace.xcm=17;
valtozok_sampletrace.xcm_blowup=5;
valtozok_sampletrace.ycm=2;
valtozok_sampletrace.ycm_current=.5;
valtozok_sampletrace.voltagelinewidth=.5;
valtozok_sampletrace.currentlinewidth=.5;
valtozok_sampletrace.markersize=1;
valtozok_sampletrace.fontsize=8;
valtozok_sampletrace.fonttype='Helvetica';
valtozok_sampletrace.axeswidth=.5;
valtozok_sampletrace.axis.voltage_y=true;
valtozok_sampletrace.axis.voltage_x=false;
valtozok_sampletrace.axis.voltage_blowup_y=false;
valtozok_sampletrace.axis.voltage_blowup_x=false;
valtozok_sampletrace.axis.current_y=false;
valtozok_sampletrace.axis.current_x=false;
valtozok_sampletrace.axis.freq_y=true;
valtozok_sampletrace.axis.freq_x=false;

persistent_generatefigures_plotselectedtimeinterval(xlsidx,dirs,xlsdata,valtozok_sampletrace,expname);


%% voltage dependence - PSD plot
ID='1411283rm_2_1,4_4';
xlsidx=find(strcmp({xlsdata.ID},ID));
expname=[ID,'_HUMAN_rhythmic_persistent_voltage_dependence'];
valtozok_PSD=struct;
valtozok_PSD.debugmode=0;
valtozok_PSD.zerotime=11838.5;
valtozok_PSD.xlimits=[-5 90]+valtozok_PSD.zerotime;
valtozok_PSD.xlimitstoexclude=[0,0]+valtozok_PSD.zerotime;
valtozok_PSD.xlimitsblowup=[25, 30;51,56;65, 70]+valtozok_PSD.zerotime;

valtozok_PSD.xcm=17;
valtozok_PSD.xcm_blowup=5;
valtozok_PSD.ycm=2;
valtozok_PSD.fontsize=8;
valtozok_PSD.fonttype='Helvetica';
valtozok_PSD.axeswidth=.5;
valtozok_PSD.caxvals=[0 .5];
valtozok_PSD.freqticks=[0.3,1,3,9];
% valtozok_PSD.freqticks=[0.5,1,2,4];

valtozok_PSD.PSDparameters=struct;
valtozok_PSD.PSDparameters.min=.25; %minimal frequency for decomposition
valtozok_PSD.PSDparameters.max=10;% - maximal frequency for decomposition
valtozok_PSD.PSDparameters.step=.01; %- frequency stepsize
valtozok_PSD.PSDparameters.scale=2;% - scale type, can be linear (1) or logarithmic (2)
valtozok_PSD.PSDparameters.wavenumber=9;% - number of waves in wavelet
valtozok_PSD.PSDparameters.waveletlength=30;
valtozok_PSD.PSDparameters.addtaper=1;
valtozok_PSD.PSDparameters.taperlength=30;

persistent_generatefigures_plot_PSD_for_timeinterval(xlsidx,dirs,xlsdata,valtozok_PSD,expname);

%%

%% stim during rhythmic PF - sample trace

ID='1612012rm_5_1,2_4';
xlsidx=find(strcmp({xlsdata.ID},ID));
expname=[ID,'_HUMAN_rhythmic_persistent_stim_during_PF'];
valtozok_sampletrace=struct;
valtozok_sampletrace.debugmode=1;
valtozok_sampletrace.zerotime=57871;
valtozok_sampletrace.xlimits=[0 120]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.ylimits=[0 0];%
valtozok_sampletrace.xlimitstoexclude=[0,0]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.xlimitsblowup=[0, 0]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.ylimitsblowup=zeros(size(valtozok_sampletrace.xlimitsblowup));
% valtozok_sampletrace.ylimitsblowup=[-72,24;-72,24;-72,24];


% ez már egy második..
ID='1605111rm_1_1_4';
xlsidx=find(strcmp({xlsdata.ID},ID));
expname=[ID,'_HUMAN_rhythmic_persistent_stim_during_PF'];
valtozok_sampletrace=struct;
valtozok_sampletrace.debugmode=1;
valtozok_sampletrace.zerotime=62100;
valtozok_sampletrace.xlimits=[0 200]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.ylimits=[0 0];%
valtozok_sampletrace.xlimitstoexclude=[0,0]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.xlimitsblowup=[0, 0]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.ylimitsblowup=zeros(size(valtozok_sampletrace.xlimitsblowup));

% ez már egy harmadik..
ID='1605102rm_4_1_4';
xlsidx=find(strcmp({xlsdata.ID},ID));
expname=[ID,'_HUMAN_rhythmic_persistent_stim_during_PF'];
valtozok_sampletrace=struct;
valtozok_sampletrace.debugmode=1;
valtozok_sampletrace.zerotime=57575;
valtozok_sampletrace.xlimits=[0 15]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.ylimits=[0 0];%
valtozok_sampletrace.xlimitstoexclude=[0,0]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.xlimitsblowup=[0, 0]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.ylimitsblowup=zeros(size(valtozok_sampletrace.xlimitsblowup));


% ez már egy negyedik..
ID='1603102rm_4_1_4';
xlsidx=find(strcmp({xlsdata.ID},ID));
expname=[ID,'_HUMAN_rhythmic_persistent_stim_during_PF'];
valtozok_sampletrace=struct;
valtozok_sampletrace.debugmode=0;
valtozok_sampletrace.zerotime=79029.7;
valtozok_sampletrace.xlimits=[0 15]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.ylimits=[0 0];%
valtozok_sampletrace.xlimitstoexclude=[0,0]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.xlimitsblowup=[2.6,3.4;9.3, 10.1]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.ylimitsblowup=zeros(size(valtozok_sampletrace.xlimitsblowup));




% valtozok_sampletrace.ylimitscurr=[0 0];%[-150,250];
valtozok_sampletrace.isiYlimits=[0 2];
valtozok_sampletrace.freqYlimits=[.1 500];
valtozok_sampletrace.freqYscale=[1, 10, 100];
valtozok_sampletrace.cutofffreq=15000;
valtozok_sampletrace.highlightaxonalspikes=1;
valtozok_sampletrace.highlightaxonalspikes_timeback=.001;
valtozok_sampletrace.highlightaxonalspikes_timeforward=.020;


valtozok_sampletrace.xcm=17;
valtozok_sampletrace.xcm_blowup=8;
valtozok_sampletrace.ycm=2;
valtozok_sampletrace.ycm_current=.5;
valtozok_sampletrace.voltagelinewidth=.5;
valtozok_sampletrace.currentlinewidth=.5;
valtozok_sampletrace.markersize=1;
valtozok_sampletrace.fontsize=8;
valtozok_sampletrace.fonttype='Helvetica';
valtozok_sampletrace.axeswidth=.5;
valtozok_sampletrace.axis.voltage_y=true;
valtozok_sampletrace.axis.voltage_x=true;
valtozok_sampletrace.axis.voltage_blowup_y=true;
valtozok_sampletrace.axis.voltage_blowup_x=false;
valtozok_sampletrace.axis.current_y=false;
valtozok_sampletrace.axis.current_x=false;
valtozok_sampletrace.axis.freq_y=true;
valtozok_sampletrace.axis.freq_x=true;


% valtozok_sampletrace.threshold.threshold=10;
% valtozok_sampletrace.threshold.stepbacklength=.002;
% valtozok_sampletrace.threshold.aplength=.006;
valtozok_sampletrace.threshold.displayonvoltagetrace=1;
% valtozok_sampletrace.threshold.displayonblowup=0;
% valtozok_sampletrace.threshold.movingstep=4;
% valtozok_sampletrace.threshold.threshmedianwindow=.15;
% valtozok_sampletrace.threshold.addbaselineval=1;

persistent_generatefigures_plotselectedtimeinterval(xlsidx,dirs,xlsdata,valtozok_sampletrace,expname);

%% pharmacology - cadmium - sample trace

ID='1605062rm_4_1_4';
xlsidx=find(strcmp({xlsdata.ID},ID));
expname=[ID,'_HUMAN_rhythmic_persistent_pharmacology_cadmium'];
valtozok_sampletrace=struct;
valtozok_sampletrace.debugmode=0;
valtozok_sampletrace.zerotime=64485;
valtozok_sampletrace.xlimits=[-80 600]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.ylimits=[-80 30];%
valtozok_sampletrace.xlimitstoexclude=[0,0]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.xlimitsblowup=[-71.8, -31.8;183.8,218.8]+valtozok_sampletrace.zerotime;
valtozok_sampletrace.ylimitsblowup=zeros(size(valtozok_sampletrace.xlimitsblowup));
% valtozok_sampletrace.ylimitsblowup=[-72,24;-72,24;-72,24];

% valtozok_sampletrace.ylimitscurr=[0 0];%[-150,250];
valtozok_sampletrace.isiYlimits=[0 2];
valtozok_sampletrace.freqYlimits=[.1 500];
valtozok_sampletrace.freqYscale=[1, 10, 100];
valtozok_sampletrace.cutofffreq=10500;
valtozok_sampletrace.highlightaxonalspikes=1;
valtozok_sampletrace.highlightaxonalspikes_timeback=.001;
valtozok_sampletrace.highlightaxonalspikes_timeforward=.020;
% valtozok_sampletrace.drugwashin.drugwashinline_ystart=30;
% valtozok_sampletrace.drugwashin.drugwashinline_ystep=5;
% valtozok_sampletrace.drugwashin.drugwashinline_linewidth=2;

valtozok_sampletrace.xcm=17;
valtozok_sampletrace.xcm_blowup=8;
valtozok_sampletrace.ycm=2;
valtozok_sampletrace.ycm_current=.5;
valtozok_sampletrace.voltagelinewidth=.5;
valtozok_sampletrace.currentlinewidth=.5;
valtozok_sampletrace.markersize=1;
valtozok_sampletrace.fontsize=8;
valtozok_sampletrace.fonttype='Helvetica';
valtozok_sampletrace.axeswidth=.5;
valtozok_sampletrace.axis.voltage_y=true;
valtozok_sampletrace.axis.voltage_x=true;
valtozok_sampletrace.axis.voltage_blowup_y=true;
valtozok_sampletrace.axis.voltage_blowup_x=true;
valtozok_sampletrace.axis.current_blowup_y=true;
valtozok_sampletrace.axis.current_blowup_x=true;
valtozok_sampletrace.axis.current_y=true;
valtozok_sampletrace.axis.current_x=true;
valtozok_sampletrace.axis.freq_y=true;
valtozok_sampletrace.axis.freq_x=true;

persistent_generatefigures_plotselectedtimeinterval(xlsidx,dirs,xlsdata,valtozok_sampletrace,expname);

