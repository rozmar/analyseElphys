% normalizeMatrixByGroup normalizes matrix individually with the containing
% cycle as baseline. Then, averages individually normalized matrices and
% returns that matrix.
%
% Parameters
%  - timeVector - tx1 vector, time values for the power matrix
%  - matrixArray - nx1 array of txm power matrices
%  - neighborMatrix - nx2 matrix, neighbor markers for normalization
% Return value
%  - normalizedMatrix - txm power matrix: normalized and averaged form of
%    the input
function normalizedMatrix = normalizeMatrixByGroup(timeVector, matrixArray, neighborMatrix, normType)

  %% -----------------------------
  %  Initialize temporary output matrix
  %% -----------------------------
  overallPowerMatrix = zeros(size(matrixArray{1}));
  %% -----------------------------
  
  %% -----------------------------
  % Normalize each matrix individually
  %% -----------------------------
  for t = 1 : length(matrixArray)
      
    % Find the interval between neighbors
    neighborInterval = (neighborMatrix(t,1)<=timeVector&...
        timeVector<=neighborMatrix(t,2));
    
    % Get the baseline power
    powerMatrix      = matrixArray{t};
    baslinePower     = powerMatrix(:,neighborInterval);
    
    % Normalize with respect to the neighbors
    powerMatrix = normalizeMatrixWithBaseline(powerMatrix, baslinePower, normType); 
    
    % Convert back from decibel scale for average
    if strcmpi(normType,'decibel')
      powerMatrix = 10.^(powerMatrix./10);    
    end
    
    % Sum powers across trials
    overallPowerMatrix = overallPowerMatrix + powerMatrix;
  end
  %% -----------------------------
 
  %% -----------------------------
  % Average matrices
  %% -----------------------------
  normalizedMatrix = overallPowerMatrix ./ t;
  %% -----------------------------
  
  %% -----------------------------
  % Convert to decibel
  %% -----------------------------
  if strcmpi(normType,'decibel')
    normalizedMatrix = 10.*log10(normalizedMatrix);
  end
  %% -----------------------------
  
end