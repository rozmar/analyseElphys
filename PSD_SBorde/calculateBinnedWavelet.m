function [resultStructure, parameters] = calculateBinnedWavelet(inputStructure, parameters)

  %% -------------------------
  %  Constant definitions
  %% -------------------------
  WAVELET = 1;
  STFFT = 2;
  %% -------------------------

  %% -------------------------
  %  Preprocess input signal
  %% -------------------------
  analysisParameters = parameters.analysis;
  inputStructure = doPreprocessing(inputStructure, analysisParameters);
  inputStructure = doResampling(inputStructure, parameters);
  
  signalStructure = getSignalStructure(inputStructure, parameters);
  
  signalStructure = doFiltering(signalStructure, parameters);
  %% -------------------------
 
  %% -------------------------
  %  Collect segment markers
  %% -------------------------
  segmentStructure = collectSegments(inputStructure, parameters);
  %% -------------------------
 
  %% -------------------------
  %  Add buffer zone to avoid edge artifact
  %% -------------------------
  bufferedSegment = addBufferZone(segmentStructure, analysisParameters);
  binnedTime = createBinnedTime(analysisParameters);
  %% -------------------------
  
  %% -------------------------
  %  Cut segments from signal
  %% -------------------------
  segmentedSignalArray = cutSegments(signalStructure, bufferedSegment);
  %% -------------------------
   
  %% -------------------------
  % Iterate segments
  %% -------------------------
  numberOfSegments = length(bufferedSegment);
  rawPowerPerSegment = cell(numberOfSegments, 1);
  rawTimeVectors = cell(numberOfSegments, 1);
  binnedPowerPerSegment = cell(numberOfSegments, 1);
  for s = 1 : numberOfSegments
    
    % Handle artefacts, if needed
    currentSignal = handleArtifacts(segmentedSignalArray{s}, analysisParameters);
    
    % This signal was discarded
    if isempty(currentSignal) 
      continue;
    end

    currentOriginalSegment = segmentStructure(s);
    currentBufferedSegment = bufferedSegment(s);
    
    % Normalize time to duration of segment    
    [currentSignalTime, segmentDuration, paddingBefore] = ...
      createTimeForSegment(currentBufferedSegment, currentOriginalSegment, length(currentSignal));
    segmentStructure(s).duration = segmentDuration;
    rawTimeVectors{s} = currentSignalTime * segmentDuration;
    
    if analysisParameters.tfdMode == WAVELET
      % Perform Wavelet TFD
      waveletParameters = analysisParameters.wavelet;
      waveletParameters.interval = getSI(signalStructure);
      
      [powerMatrix, frequencyVector] = ...
       calculateTFDPower(currentSignal, waveletParameters);
      timeVector = currentSignalTime;
    
    elseif analysisParameters.tfdMode == STFFT
      % Perform Short-time Fourier transform
      stftParameters = analysisParameters.stfft;
      
      Fs = getFS(signalStructure);
      window = round(stftParameters.windowSize * Fs);
      overlap = round((stftParameters.windowSize - stftParameters.shiftSize) * Fs);
      freqs = createStfftScale(stftParameters);
      [FFT, frequencyVector, timeVector] = spectrogram(currentSignal, window, overlap, freqs, Fs);
      
      powerMatrix = abs(FFT) .^ 2;  % Calculate power from coefficient
      timeVector = (toCol(timeVector) - paddingBefore) ./ segmentDuration;
    
    end

    % Save power matrix without binning 
    rawPowerPerSegment{s} = powerMatrix;
    
    % Bin and average power values
    binnedPowerPerSegment{s} = averageAndBinMatrixByTime(powerMatrix, timeVector, binnedTime);
    
    if isSinglePacketPlot(parameters.plot)
      % Plot the segment trace and spectrogram for manual validation
      doPlotWithSignal(timeVector, frequencyVector, powerMatrix, currentSignal, currentSignalTime, parameters);
      title(sprintf('Segment between %f - %f', currentOriginalSegment.start, currentOriginalSegment.end));
      subplot(2,1,1);
      title(sprintf('Raw segment, %.3f s', currentOriginalSegment.end - currentOriginalSegment.start));
    end
  end
  %% ---------------------
  
  if isAlignedAveragePlot(parameters.plot)
    alignedStructure = alignAndAveragePowers(segmentStructure, rawPowerPerSegment, rawTimeVectors, frequencyVector, parameters.plot);
    subplot(2, 1, 1);
    title(sprintf('Averaged aligned wavelet for %s', unescape(inputStructure.title)));
  end  
  
  %% -------------------------
  %  Remove discarded segments
  %% -------------------------
  discardedSegment = cellfun(@isempty, binnedPowerPerSegment);
  bufferedSegment(discardedSegment) = [];
  binnedPowerPerSegment(discardedSegment) = [];
  numberOfSegments = length(binnedPowerPerSegment);
  
  % If nothing remained, create an empty return value
  if numberOfSegments == 0
    resultStructure.powerPerSegment = {};
    resultStructure.averagedPower = [];
    resultStructure.frequencyVector = [];
    resultStructure.segmentStructure = [];
    resultStructure.timeVector = [];
    return;
  end
  %% -------------------------
  
  %% ---------------------
  %  Average the power matrices 
  %% ---------------------
  averagedPowerMatrix = nanmean_matrices(binnedPowerPerSegment);
  
  % Remove last value, it is only for binning
  binnedTime(end) = [];
  %% ---------------------

  %% ---------------------
  %  Plot average if needed
  %% ---------------------
  if isSingleCellAveragePlot(parameters.plot)
    doPlotWithSignal(binnedTime, frequencyVector, averagedPowerMatrix, [], [], parameters);
    title(sprintf('Averaged binned wavelet for %s', unescape(inputStructure.title)));
  end
  %% ---------------------
 
  %% ---------------------
  % Store the results
  %% ---------------------
  resultStructure.powerPerSegment = binnedPowerPerSegment;
  resultStructure.averagedPower = averagedPowerMatrix;
  resultStructure.frequencyVector = frequencyVector;
  resultStructure.segmentStructure = bufferedSegment;
  resultStructure.timeVector = binnedTime;  
  if exist('alignedStructure', 'var')
    resultStructure.aligned = alignedStructure;
  end
  %% ---------------------
  
end

%% ========================
%  Functions for binning
%% ========================
% Create normalized binned time
function binnedTime = createBinnedTime(parameters)

  % Size of buffer (or ratio, depends on mode)
  bufferSize = parameters.bufferSize;
  % Number of bins during segment
  numOfBins = parameters.numberOfBins;

  % Calculate the number of bins at padding
  if parameters.bufferMode == 2
    paddingBinNum = round(bufferSize * numOfBins);
  else
    paddingBinNum = round(numOfBins * bufferSize / parameters.minSegmentLength);
    bufferSize = bufferSize / parameters.minSegmentLength;
  end
  
  % Create binned time vector for binnnig.
  binnedTime = unique([ ...
    linspace(-bufferSize, 0, paddingBinNum + 1), ...
    linspace(0, 1, numOfBins + 1), ...
    linspace(1, 1+bufferSize, paddingBinNum + 1)])';
end
%% ========================

%% ========================
%  Buffer zone functions
%% ========================
% Add buffer to each segment
function segmentStructure = addBufferZone(segmentStructure, parameters)

  bufferSize = parameters.bufferSize; 
  bufferMode = parameters.bufferMode;

  % Put parameter to correct variable
  if bufferMode == 1
    bufferLength = bufferSize;
  elseif bufferMode == 2
    bufferRatio = bufferSize;
  end

  % Iterate segments and add buffer to them
  numberOfSegment = length(segmentStructure);
  for i = 1 : numberOfSegment
    
    if bufferMode == 2
      bufferLength = (segmentStructure(i).end - segmentStructure(i).start) * bufferRatio;  
    end
    
    segmentStructure(i) = addBufferToSegment(segmentStructure(i), bufferLength);
  end
  
end

% Add buffer zone of given length to a given segment 
function segment = addBufferToSegment(segment, bufferLength)
  segment.start = max(segment.start - bufferLength, 0);
  segment.end = segment.end + bufferLength;
end
%% ========================

%% ========================
%  Segment handling functions
%% ========================
% Collects segments from input and returns in the new format
% Eliminates those segment which were filtered out by the given predicates
function segmentStructure = collectSegments(inputStructure, parameters)

  segmentName = parameters.channel.segment;

  if strcmpi(segmentName.name, 'whole')
    segmentStructure.start = 0;
    segmentStructure.end = inputStructure.(parameters.channel.signal).times(end);
    return;
  end

  segmentStructure = getSingleSegment(inputStructure, segmentName);
  segmentStructure = convertSegmentNewFormat(segmentStructure);
  segmentStructure = doSegmentElimination(segmentStructure, parameters.analysis);
end

% Create a time vector for a segment with a given number of samplepoint
function [timeVector, segmentDuration, paddingBefore] = createTimeForSegment(bufferedSegment, originalSegment, samplePoint)
  % Calculate original segment duration
  segmentDuration = originalSegment.end - originalSegment.start;
  % Generate time vector for buffered segment
  timeVector = linspace(bufferedSegment.start, bufferedSegment.end, samplePoint)';
  % Normalize timevector according to the original segment
  timeVector = (timeVector - originalSegment.start) ./ segmentDuration;
  % Original padding value
  paddingBefore = originalSegment.start - bufferedSegment.start;
end
%% ========================

%% ========================
%  Utility functions
%% ========================
% Handle artefact in signal. Handling method can be simple discard or
% outlier removal. Different outlier removal methods can be applied. If
% discard mode have been selected, the whole signal will be discarded if
% outlier value exists. If no handling mode was given, artefacts will be
% retained.
function signal = handleArtifacts(signal, parameters)

  if ~isfield(parameters, 'artefactHandling')
    return;
  end

  if strcmp(parameters.artefactHandling, 'discard')
    
    outlierValues = detectOutliers(signal, 'abs', 5);
    
    if ~isempty(outlierValues)
      signal = [];
    end
    
  elseif strcmp(parameters.artefactHandling, 'remove')
    signal = removePotentialArtefact(signal, parameters);
  end

end

% Decides whether we need to plot the single segment plot
function answer = isSinglePacketPlot(parameters)
  answer = isfield(parameters, 'singlePacket') && parameters.singlePacket;
end

% Decides whether we need to plot the aligned segment average plot
function answer = isAlignedAveragePlot(parameters)
  answer = isfield(parameters, 'singleCellAlignedAverage') && parameters.singleCellAlignedAverage;
end

% Get the trigger of alignment of the aligned average plot. 
function answer = getAlignedAveragePlotTrigger(parameters)
  if isfield(parameters, 'singleCellAlignedAverageTrigger')
    answer = parameters.singleCellAlignedAverageTrigger;
  else
    answer = 'start';
  end
end

% Decides whether we need to plot the average segment plot
function answer = isSingleCellAveragePlot(parameters)
  answer = isfield(parameters, 'singleCellAverage') && parameters.singleCellAverage;
end

% Unescape filename for plot title
function unescaped = unescape(originalString)
  unescaped = strrep(originalString, '_', '');
end

% Returns the signal structure from the input by its name
function signalStructure = getSignalStructure(inputStructure, parameters)
  signalName = parameters.channel.signal;
  signalStructure = inputStructure.(signalName);
end

% Create frequency vector with given scaling (linear or log)
function freqs = createStfftScale(stftParameters)
  if stftParameters.scale == 1
    freqs = linspace(stftParameters.min, stftParameters.max, stftParameters.numFreq+1);
  elseif stftParameters.scale == 2
    freqs = logspace(log10(stftParameters.min), log10(stftParameters.max), stftParameters.numFreq+1);
  end
end

% Do spectrogram plot
function doPlotWithSignal(timeVector, frequencyVector, powerMatrix, signalVector, signalTime, parameters)

  withSignal = ~isempty(signalVector);

  figure;
  if withSignal
    subplot(2,1,1);
  end
  
  if parameters.analysis.tfdMode == 1 && isfield(parameters.plot, 'timeBound')
    tB = parameters.plot.timeBound;
    remTime = (tB(1)<=timeVector&timeVector<=tB(2));
    timeVector(~remTime) = [];
    powerMatrix(:, ~remTime) = [];
    %powerMatrix = nanzscore(powerMatrix);
  else
    %powerMatrix = log10(powerMatrix);
  end
  
  plotPowerSpectrum(timeVector, frequencyVector, powerMatrix, parameters.plot);
  %plotPowerSpectrum(timeVector, frequencyVector, zscore(log10(powerMatrix)')', parameters.plot);
  %plotPowerSpectrum(timeVector, frequencyVector, nanzscore(powerMatrix), parameters.plot);
  %plotPowerSpectrum(timeVector, frequencyVector, powerMatrix, parameters.plot);
  %title(sprintf('Raw segment, %.3f s', timeVector(end)));
  
  if withSignal
    xLimit = get(gca, 'xlim');
  
    subplot(2,1,2);
    plot(signalTime, signalVector);
    xlim(xLimit);
  end
  
end
%% ========================
