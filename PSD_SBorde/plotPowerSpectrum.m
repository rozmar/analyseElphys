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
  %  Make restrictions before plotting
  %% ----------------------------
  [timeVector, powerMatrix, xLimit] = restrictByTime(timeVector, powerMatrix, parameters);
  [frequencyVector, powerMatrix, yLimit] = restrictByFrequency(frequencyVector, powerMatrix, parameters);
  %% ----------------------------   
  
  %% ----------------------------
  %  Plot spectrogram
  %% ----------------------------  
  pcolor(timeVector, frequencyVector, powerMatrix);
  colormap jet;
  shading interp;
  %% ----------------------------
  
  %% ----------------------------
  %  Customize plot
  %% ----------------------------  
 
  % Set axis limit
  setXAxis(xLimit, parameters);
  setYAxis(yLimit, parameters);
  setZAxis(parameters);
  
  % Set title for plot
  createTitle(parameters);
  %% ----------------------------  

end

%% ---------------------------------
%  Axis customization function
%% ---------------------------------
% Customize the X axis
function setXAxis(xLimit, parameters)

  if isfield(parameters, 'xTick')
    set(gca, 'xtick', parameters.xTick);
    set(gca, 'xticklabel', parameters.xTickLabel);
  end
  
  set(gca, 'xlim', xLimit);
end

% Customize the Y axis
function setYAxis(yLimit, parameters)

  % Create ticks
  if ~isfield(parameters, 'yTick') || isempty(parameters.yTick)
    if isLogScale(parameters)
      yticks = logspace(log10(yLimit(1)), log10(yLimit(2)), 10);
    elseif isLinScale(parameters)
      yticks = linspace(yLimit(1), yLimit(end), 10);
    end
  else
    yticks = parameters.yTick; 
  end
  
  % Set scale type
  if isLogScale(parameters)
    set(gca, 'yscale', 'log');
  elseif isLinScale(parameters)
    set(gca,'yscale','lin');
  end
  
  % Add ticks, ticklabels and set limit
  set(gca, 'ytick', yticks);
  set(gca, 'yticklabel', round(yticks*10)/10); 
  set(gca, 'ylim', yLimit);
end

% Customize Z axis (i.e. the color map)
function setZAxis(parameters)
  
  if isfield(parameters, 'colorLimit') && ~isempty(parameters.colorLimit)
    set(gca, 'clim', parameters.colorLimit);
  end
  
end
%% ---------------------------------

%% ---------------------------------
%  Value restriction before plot
%% ---------------------------------
% This function performs the restriction in time (if restriction needed)
function [timeVector, powerMatrix, xLimit] = restrictByTime(timeVector, powerMatrix, parameters)

  timeVector = createTimeVector(timeVector, size(powerMatrix, 2));

  if isfield(parameters, 'timeBound')
    powerMatrix = restrictValue(powerMatrix, timeVector, parameters.timeBound);
    timeVector = restrictValue(toRow(timeVector), timeVector, parameters.timeBound);
  end
  
  xLimit = timeVector([1,end]);
end

% Creates time vector if it wasn't given, else returns that
function timeVector = createTimeVector(timeVector, numberOfSamples)
  if isempty(timeVector)
    timeVector = (1:numberOfSamples)'; 
  end
end

% This function performs the restriction in frequency (if restriction needed)
function [frequencyVector, powerMatrix, yLimit] = restrictByFrequency(frequencyVector, powerMatrix, parameters)
  
  if isfield(parameters, 'frequencyBound')
    powerMatrix = restrictValue(powerMatrix', frequencyVector, parameters.frequencyBound)';
    frequencyVector = toCol(restrictValue(toRow(frequencyVector), frequencyVector, parameters.frequencyBound));
  end
  
  yLimit = frequencyVector([1,end]);
  
end

% Find the values in restrictionValue vector between the permitted limit,
% and remove elements outside this limit from the value matrix.
function restrictedValue = restrictValue(originalValue, restrictionValue, limit)
  indices = findRemainingIndices(restrictionValue, limit);
  restrictedValue = originalValue(:, indices);
end

% Returns the indices of the values from a vector which falls between the
% given limits.
function indices = findRemainingIndices(value, limit)
  indices = (limit(1) <= value & value <= limit(2));
end
%% ---------------------------------

%% ---------------------------------
%  Utility functions
%% ---------------------------------
% Decides whether we should use linear scale
function answer = isLinScale(parameters)
  answer = ~isfield(parameters, 'scale') || parameters.scale == 1;
end

% Decides whether we should use log scale
function answer = isLogScale(parameters)
  answer = isfield(parameters, 'scale') && parameters.scale == 2;
end

% Set plot title
function createTitle(parameters)
  if isfield(parameters, 'title')
    titleString = parameters.title;
  else
    titleString = 'Spectrogram';
  end
  
  title(titleString);
end
%% ---------------------------------
