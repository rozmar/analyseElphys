function [dirs,xlsdata,projectdata,amplifier]=aE_projectselector_main
% aE_projectselector_main generates variables for analysElphys_main.m
% this script is instrumental for proper functiong of analysElphys_main.m
% edit this file if you would like to add further projects (expand both the 
% projectnames variable and the if-elseif structure, projectnum refers to the 
% index of the selected projectname in the projectnames list)
% the outputs are:
%       - dirs - the directory structure of the selected project
% with the location of raw data, bridged data, events, power spectrum
% densities etc.
%       - xlsdata - the metadata entered in an .xls file by the
% user
%       - projectdata - metadata entered by the user thorugh a 
% GUI (aE_projectselector.m). this contains information about whichs scripts
% are to run in the analysElphys_main.m script
%       - amplifier - HEKA or axon
projectnames={'CB1elphys','InVivo','Persistent-ChRstim','persistent firing','bleb recording'}; % 
projectdata.owbridge=0;
projectdata.owbridge=0;
projectdata.owevent=0;
projectdata = aE_projectselector(projectnames);
projectnum=projectdata.projectnum;
% alapadatok
if projectnum==1; % rosehip CB1 experiments
    overwrite=0;
    locations=marcicucca_locations;
    dirs.basedir=[locations.tgtardir,'ANALYSISdata/marci/Human_rosehip/CB1elphys/'];
    dirs.rawexporteddir=[dirs.basedir,'Exported_raw/'];
    dirs.bridgeddir=[dirs.basedir,'Bridged_stim/'];
    dirs.eventdir=[dirs.basedir,'Events/'];
    dirs.eventparaleldir=[dirs.basedir,'Events/paralel/'];
    dirs.onlyAPeventdir=[dirs.basedir,'Events_onlyAP/'];
    dirs.grpupedeventdir=[dirs.basedir,'Events_grouped/'];
    dirs.stimepochdir=[dirs.basedir,'Stimepochs/'];
    dirs.figuresdir=[dirs.basedir,'Figures/'];
    amplifier='HEKA';
    xlsdata=aE_readxls([dirs.basedir,'cb1elphys.xls']);
elseif projectnum==2; % in vivo layer 1
    overwrite=0;
    locations=marcicucca_locations;
    dirs.basedir=[locations.tgtardir,'ANALYSISdata/marci/_persistent/_InVivo/'];
    %     dirs.rawdir=[locations.tgtardir,'AXONdata/'];
    dirs.bridgeddir=[dirs.basedir,'Bridged_stim/'];
    dirs.rawexporteddir=[dirs.basedir,'Exported_raw/'];
    dirs.eventdir=[dirs.basedir,'Events/'];
    dirs.eventparaleldir=[dirs.basedir,'Events/paralel/'];
    dirs.statedir=[dirs.basedir,'State/'];
    dirs.taxonomydir=[locations.tgtardir,'ANALYSISdata/marci/_Taxonomy/persistent_invivo/'];
    dirs.PSDdir=[dirs.basedir,'PSD/'];
    dirs.PSDdir_high=[dirs.basedir,'PSD_high/'];
    dirs.PSDdir_log=[dirs.basedir,'PSD_log/'];
    dirs.videodir=[dirs.basedir,'Videodata/'];
    dirs.brainstatedir=[dirs.basedir,'BrainState/'];
    % dirs.onlyAPeventdir=[dirs.basedir,'Events_onlyAP/'];
    % dirs.grpupedeventdir=[dirs.basedir,'Events_grouped/'];
    % dirs.stimepochdir=[dirs.basedir,'Stimepochs/'];
    dirs.figuresdir=[dirs.basedir,'Figures/'];
    dirs.breathingdir=[dirs.basedir,'Breathing/'];
    dirs.breathingdir_psd=[dirs.basedir,'Breathing_PSD/'];
    dirs.offsetdir=[dirs.basedir,'Offset/'];
    dirs.cross_spectrum_breathing=[dirs.basedir,'Cross_spectrum_with_breathing/'];
    dirs.NWBdir=[dirs.basedir,'NWB/'];
    xlsdata=aE_readxls([dirs.basedir,'invivodata.xls']);
    amplifier='HEKA';
elseif projectnum==3; % ChannelRodopsin stimulaion experiments for axonalAP inducion
    overwrite=0;
    locations=marcicucca_locations;
    dirs.basedir=[locations.tgtardir,'ANALYSISdata/marci/_persistent/_ChRstim/'];
    dirs.rawexporteddir=[dirs.basedir,'Exported_raw/'];
    dirs.bridgeddir=[dirs.basedir,'Bridged_stim/'];
    dirs.eventdir=[dirs.basedir,'Events/'];
    dirs.eventparaleldir=[dirs.basedir,'Events/paralel/'];
    dirs.PSDdir=[dirs.basedir,'PSD/'];
    %     dirs.onlyAPeventdir=[dirs.basedir,'Events_onlyAP/'];
    %     dirs.grpupedeventdir=[dirs.basedir,'Events_grouped/'];
    %     dirs.stimepochdir=[dirs.basedir,'Stimepochs/'];
    dirs.figuresdir=[dirs.basedir,'Figures/'];
    amplifier='HEKA';
    xlsdata=aE_readxls([dirs.basedir,'ChRstimdata_windows.xls']);
elseif projectnum==4 % persistent firing human and rodent brain slice experiments
    overwrite=0;
    locations=marcicucca_locations;
    dirs.basedir=[locations.EMdir,'ANALYSISdata/marci/_persistent/'];
    dirs.rawexporteddir=[dirs.basedir,'Exported_raw/'];
    dirs.bridgeddir=[dirs.basedir,'Bridged_stim/'];
    dirs.taxonomydir=[locations.tgtardir,'ANALYSISdata/marci/_Taxonomy/persistent_slice/'];
    dirs.eventdir=[dirs.basedir,'Events/'];
    dirs.eventparaleldir=[dirs.basedir,'Events/paralel/'];
    dirs.onlyAPeventdir=[dirs.basedir,'Events_onlyAP/'];
    dirs.PSDdir=[dirs.basedir,'PSD/'];
    dirs.PSDdir_log=[dirs.basedir,'PSD_log/'];
    dirs.grpupedeventdir=[dirs.basedir,'Events_grouped/'];
    dirs.stimepochdir=[dirs.basedir,'Stimepochs/'];
    dirs.figuresdir=[dirs.basedir,'figures/'];
    amplifier='HEKA';
    xlsdata=aE_readxls([dirs.basedir,'persistentdata_windows.xls']);
    dirs.v0distdir=[dirs.basedir,'v0_dist/'];
    dirs.v0_vs_PSDdir=[dirs.basedir,'v0_vs_PSD/'];
elseif projectnum==5 % cell attached bleb recordings
    overwrite=0;
    locations=marcicucca_locations;
    dirs.basedir=[locations.tgtardir,'ANALYSISdata/marci/_persistent/_BlebRecording/'];
    dirs.rawexporteddir=[dirs.basedir,'Exported_raw/'];
    dirs.bridgeddir=[dirs.basedir,'Bridged_stim/'];
    dirs.eventdir=[dirs.basedir,'Events/'];
    dirs.eventparaleldir=[dirs.basedir,'Events/paralel/'];
    %     dirs.onlyAPeventdir=[dirs.basedir,'Events_onlyAP/'];
    %     dirs.grpupedeventdir=[dirs.basedir,'Events_grouped/'];
    %     dirs.stimepochdir=[dirs.basedir,'Stimepochs/'];
    dirs.figuresdir=[dirs.basedir,'Figures/'];
    amplifier='HEKA';
    xlsdata=aE_readxls([dirs.basedir,'blebdata_windows.xls']);
end