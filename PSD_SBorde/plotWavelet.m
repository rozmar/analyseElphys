function plotWavelet(resultStructure, parameters)
  %% ------------------------------
  %  Normalize and plot power spect.
  %% ------------------------------ 
  if parameters.plot.overallAverage
      
    if isGlobalBaseline(parameters)
      normPower = normalizeSpectrogramGlobal(powerMatrix, displayInterval, baselineInterval, parameters);  
    else
      normPower = normalizeSpectrogramIndividual(pPerTrial, displayInterval, timeVector, neighborMatrix, parameters);
    end
    
    trText   = parameters.channel.triggerTitle;
    normType = parameters.wavelet.normalization.type;
    
    parameters.plot.title = sprintf('%s triggered %s (%s norm.)', trText, strrep(fileName,'_','\_'), normType);
    
    figure;
    displaySpectrogram(timeVector(displayInterval), frequencyVector, normPower, sliceMatrix(:,displayInterval), parameters.plot);
  end
  
  if parameters.plot.singlePlot
      
    nTrial = length(pPerTrial);
    
    for t = 1 : nTrial
      normalizeAndDisplaySpectrogram(pPerTrial{t}, sliceMatrix(t,:), displayInterval, baselineInterval, timeVector(displayInterval), frequencyVector, fileName, parameters);
    end
  end
end