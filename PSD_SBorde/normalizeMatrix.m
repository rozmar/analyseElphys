% normalizeMatrix normalize each row of a given matrix. Normalization
% method can be given in parameter. (Default is the z-score).
%
% normalizedMatrix = normalizeMatrix(inputMatrix, normalizationType)
%   inputMatrix - a matrix to normalize
%   normalizationType - type of normalization (optional)
%                        1 - z-score (default)
%                        2 - decibel
%   normalizedMatrix - normalized matrix
function normalizedMatrix = normalizeMatrix(inputMatrix, normalizationType)

  if nargin<2
    normalizationType = 1;
  end
  
  if strcmpi(normalizationType,'zscore')
    % Z-score normalization
    normalizedMatrix = zscore(inputMatrix')';
  elseif strcmpi(normalizationType,'decibel')
    % Decibel normalization
    meanMatrix = mean(inputMatrix, 2)*ones(1, size(inputMatrix,2));
    normalizedMatrix = 10.*log10(inputMatrix ./ meanMatrix);
  end
end