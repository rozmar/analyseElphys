% normalizeMatrixWithBaseline normalizes each row of a given matrix.
% Normalization method can be specified in parameter. The default is the
% z-score normalization. The normalization will be performed relative to a
% given baseline matrix.
%
% Parameters
%  - inputMatrix - nxt matrix, raw power have to be normalized 
%  - baseLineMatrix - nxb matrix with the same row number as inputMatrix,
%    this will be the reference of the normalization
%  - normalizationType - type of the normalization (optional), value can be
%     - 'zscore' - (default), subtract the mean of the baseline from the
%     raw matrix, then divides it by the SD of baseline power
%     - 'decibel' - divides the raw matrix by the mean of the baseline,
%     then take it's 10th logarithm and multiply by 10
% Return value
%  - normalizedMatrix - nxt matrix, the normalized form of the input
function normalizedMatrix = normalizeMatrixWithBaseline(inputMatrix, baseLineMatrix, normalizationType)

  %% ------------------------------
  %  Set default value
  %% ------------------------------
  if nargin<3
    normalizationType = 'zscore';
  end
  %% ------------------------------
  
  %% ------------------------------
  %  Do the normalization
  %% ------------------------------
  meanByRow  = mean(baseLineMatrix,2);
  colNumber  = size(inputMatrix,2);
  meanMatrix = meanByRow * ones(1, colNumber);
  if strcmpi(normalizationType,'zscore') % Z-score normalization
    sdByRow   = std(baseLineMatrix, 0, 2);
    stdMatrix = sdByRow * ones(1, colNumber);
    normalizedMatrix = (inputMatrix - meanMatrix)./stdMatrix;
  elseif strcmpi(normalizationType,'decibel') % Decibel normalization
    normalizedMatrix = 10.*log10(inputMatrix ./ meanMatrix);
  end
  %% ------------------------------
end