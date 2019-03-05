% plotWaveform plots a waveform on spectrogram. In parameter, expects the
% matrix of waveforms, and plot their average on the current figure.
% Scaling and positioning will be performed according to the parameters.
% 
% Parameters
%  - timeVector - tx1 vector, the x-axis of the plot
%  - sliceMatrix - kxt matrix, waveform of k trial
%  - parameters - customization parameters with the following fields
%     - frequencyBound - 1x2 vector, ends of the interval
%  - plotRange  - 1x2, ends of the interval between we want to strech the
%    waveform
function plotWaveform(timeVector, sliceMatrix, parameters, plotRange)

  %% ------------------------
  %  Calculate mean waveform
  %% ------------------------
  meanWaveform = mean(sliceMatrix,1);               %average trials
  meanWaveform = meanWaveform-min(meanWaveform);    %pull to zero
  waveYRange   = max(meanWaveform);                 %calculate range
  %% ------------------------
  
  %% ------------------------
  %  Shift waveform 
  %% ------------------------  
  if nargin<4
    plotYRange = diff(parameters.frequencyBound)/2;
    offset     = plotYRange/2;
  else
    plotYRange = diff(plotRange);
    offset     = plotRange(1);
  end
  
  multiplier = plotYRange/waveYRange;
  %% ------------------------ 
  
  %% ------------------------ 
  %  Plot waveform
  %% ------------------------ 
  hold on;
  plot(timeVector, multiplier*meanWaveform+offset, 'w-', 'linewidth', 2);
  hold off;
  %% ------------------------ 
  
end