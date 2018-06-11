% calculateResampledAverage calculates a resampled averaged matrix from more
% matrices, based on a given common time vector. First, resamples each
% matrix to the samplingrate of the given common time vector.
% After resample averages matrices.
% 
% Parameters
%  - matrixStructure - kx1 structure array which contains the matrices 
%    which have to be resampled and averaged. The structure have to contain 
%    the following fields:
%     - timeVector - n_1x1 vector, the time axis of the whole interval
%     - displayInterval - n_1x1 logical vector, indicates the timepoint in
%       the display interval
%     - normPower   - kxn_1_1 matrix, the power matrices to resample and
%       average, it contains only the values to display
%  - commonTime - a mx1 vector, which contains the new time vector on which
%    the resampling will be based
% Return values
%  - averageMatrix - kxm matrix, the average of the resampled matrices
function averageMatrix = calculateResampledAverage(matrixStructure, commonTime)

  %% ------------------------------
  %  Initialization
  %% ------------------------------
  nMatrix       = length(matrixStructure);
  nFrequency    = size(matrixStructure(1).normPower, 1);
  lengthTime    = length(commonTime);
  averageMatrix = zeros(nFrequency, lengthTime);    %output matrix
  %% ------------------------------
  
  %% ------------------------------
  %  Take each frequency band
  %% ------------------------------
  warning('off', 'MATLAB:linearinter:noextrap');
  for i = 1 : nFrequency
      
    currentMean = zeros(1, lengthTime);
    
    %% ------------------------------
    %  Take each matrix
    %% ------------------------------
    for j = 1 : nMatrix
        
      thisStructure = matrixStructure(j);
      thisRow = thisStructure.normPower(i,:);
      if isfield(thisStructure, 'displayTimeVector')
        thisTime = thisStructure.displayTimeVector;
      else
        thisTime = thisStructure.timeVector(thisStructure.displayInterval);
      end
       
      currentLine = timeseries(thisRow', thisTime);
      currentLine = resample(currentLine, commonTime);
      currentMean = currentMean + currentLine.data';
    end
    %% ------------------------------
    
    currentMean = currentMean ./ nMatrix;
    averageMatrix(i,:) = currentMean;
  end
  warning('on', 'MATLAB:linearinter:noextrap');
  %% ------------------------------
end
