% plotPowerSpectrum displays a power spectrogram of a given power matrix.
%
% Parameters
%  - timeVector - nx1 vector contains time, x axis
%  - frequencyVector - mx1 vector contains pseudo-frequency, y axis
%  - powerMatrix - mxn matrix, contains (normalized) power to plot, z axis
%  - parameters - plot configuring parameters with the following fields. If
%  a field wasn't given, default value will be used.
%     - scale - scale of the y axis, 1 linear (default), 2 logarithmic
%     - timeBound - time constraint for the plot, default the whole vector
%     - frequencyBound - frequency constraint for the plot, default the whole vector
%     - colorLimit - limits for the coloring, default automatic scaling
%     - title - title of the plot
function plotPowerSpectrum(timeVector, frequencyVector, powerMatrix, parameters)

  %% ----------------------------
  %  Parameter settings 
  %% ----------------------------
  %  Get scale type
  %   1 - linear
  %   2 - logarithmic
  if isfield(parameters, 'scale')
    scaleType = parameters.scale;
  else
    scaleType = 1;
  end
  
  % Cut down time
  if isfield(parameters, 'timeBound')
    timeBound   = parameters.timeBound;
    timeIndices = (timeBound(1)<=timeVector&timeVector<=timeBound(2));
    timeVector  = timeVector(timeIndices);
    powerMatrix = powerMatrix(:,timeIndices);
  else
    timeBound   = timeVector([1,end]);
  end
  
  % Cut down frequency
  if isfield(parameters, 'frequencyBound')
    freqBound        = parameters.frequencyBound;
    frequencyIndices = (freqBound(1)<=frequencyVector&frequencyVector<=freqBound(2));
    frequencyVector  = frequencyVector(frequencyIndices);
    powerMatrix      = powerMatrix(frequencyIndices, :);
  else
    freqBound        = frequencyVector([1,end]);
  end
  
  % Set title
  if isfield(parameters, 'title')
    plotTitle = parameters.title;
  else
    plotTitle = 'Spectrogram';
  end
  
  % Set color limits
  if isfield(parameters, 'colorLimit') && ~isempty(parameters.colorLimit)
    colorLimit = parameters.colorLimit;  
  else
    colorLimit = [min(powerMatrix(:)), max(powerMatrix(:))];
  end
  %% ----------------------------   
  
  %% ----------------------------
  %  Plot spectrogram
  %% ----------------------------  
  pcolor(timeVector, frequencyVector, powerMatrix);
  colormap jet;
  shading interp;
  %% ----------------------------
  
  %% ----------------------------
  %  Configure plot
  %% ----------------------------  
  if scaleType==2
    yticks = logspace(log10(frequencyVector(1)),log10(frequencyVector(end)),10);
    set(gca, 'yscale', 'log');
  elseif scaleType==1
    if ~isfield(parameters, 'yTick') || isempty(parameters.yTick)
      yticks = linspace(frequencyVector(1),frequencyVector(end),10);
    else
      yticks = parameters.yTick; 
    end
    set(gca,'yscale','lin');
  end
  %% ----------------------------
        
  set(gca, 'ytick',      yticks);
  set(gca, 'yticklabel', round(yticks*10)/10);  
  set(gca, 'clim',       colorLimit);
  
  if isfield(parameters,'xTick')
    set(gca, 'xtick',      parameters.xTick);
    set(gca, 'xticklabel', parameters.xTickLabel);
  end
  
  set(gca, 'xlim', timeBound);
  set(gca, 'ylim', freqBound);
  
  title(plotTitle);

end