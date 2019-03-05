% calculateTFDPower performs the time-frequency decomposition with complex
% Morlet wavelet with the given parameters. Returns the power matrix. It
% can handle one or more trials. If more trials were given, they will be 
% averaged. If the 4. output parameter was given, it returns the power matrix
% for every single trial too. 
%  
%
% Parameters
%  - inputSignal - txn matrix, where t is the number of samplepoints in the
%  signal, n is the number of trials. If n=1, no averaging need to be
%  performed
%  -  parameters - parameter structure with the following field
%      min - minimal frequency for decomposition
%      max - maximal frequency for decomposition
%      step - frequency stepsize
%      scale - scale type, can be linear (1) or logarithmic (2)
%      wavenumber - number of waves in wavelet
%      interval - sample interval for the signal
%      waveletlength - length of the wavelet vector in seconds - default value is 1s;
%      addtaper - (1) to add taper to avoid edge artefact;
%       taperlength - length of taper in seconds;
% Return values
%   - PowerMatrix fxt matrix, the calculated power values, f is the number 
%   of frequencies
%   - frequencyVector fx1 vector, contains the frequencies which has been
%    used
%   - CoeffMatrix fxt matrix contains the Coefficients (from which the
%    power has been calculated)
%   - PowerMatrixTrial - tx1 cell array, the power matrices for each trial
function [PowerMatrix, frequencyVector, CoeffMatrix, PowerMatrixTrial] = calculateTFDPower(inputSignal, parameters)

  %% -------------------------
  %  Constant definitions
  %% -------------------------
  LINEAR = 1;
  %% -------------------------

  %% -------------------------
  %  Set scale type (lin or log)
  %% -------------------------
  scale = LINEAR;
  
  if isfield(parameters, 'scale')
    scale = parameters.scale;
  end
  %% -------------------------

  %% -------------------------
  %  Parameter setting
  %% -------------------------  
  % Can we average the power in place?
  canAverage = ...
      isfield(parameters, 'normalization') && ...
      isfield(parameters.normalization, 'baseline') && ...
      regexpi(parameters.normalization.baseline, '(trigger|window)');
  
  % Do we need the power per trial?
  needPerTrial = (nargout==4);
 
  % If inputSignal is matrix, there is more trial
  [r,c]         = size(inputSignal);
  hasMultiTrial = (r>1&&c>1);  
  
  waitbarMsg = 'Wavelet transform %.2fHz/%.2fHz';
  
  freqRange       = [parameters.min,parameters.max];   
  frequencyVector = createScale(freqRange, parameters.step, scale);
  nFrequency      = length(frequencyVector);
  %% -------------------------   
  
  %% ------------------------- 
  %  Initialization
  %% -------------------------
  %  If we can average the power
  %  create output matrix. Else,
  %  we have to create cell array.
  %% -------------------------
  if canAverage || ~hasMultiTrial
      
    PowerMatrix = zeros(nFrequency, r);
    CoeffMatrix = zeros(nFrequency, r);
    
    if needPerTrial
      PowerMatrixTrial = cell(c,1);
    end
    
  elseif regexpi(parameters.normalization.baseline, 'interval')
      
    PowerMatrix      = [];
    CoeffMatrix      = [];
    PowerMatrixTrial = cell(c,1);
    
  end
  %% -------------------------
    
  %% -------------------------
  %  In case of more trial, 
  %  have to be reshaped
  %% -------------------------
  if hasMultiTrial
    inputSignal = reshape(inputSignal, 1, r*c);
  end
  %% -------------------------
    
  %% -------------------------
  %  Process each frequency
  %% -------------------------
  for i = 1 : nFrequency
      
    f = frequencyVector(i);
       
    % Display waitbar
    newTitle = sprintf(waitbarMsg, frequencyVector([i,end]));
    if i==1
%       multiWaitbar( newTitle , 'Value',  i/nFrequency);  
    elseif i > 1
        %       multiWaitbar( lastTitle , 'Value', i/nFrequency, 'Relabel', newTitle);
    end
    lastTitle = newTitle;
    
    % Create wavelet and perform convolution
    wavelet  = createMorletWavelet(f, parameters);
    if isfield(parameters,'addtaper')
        taperlength=min(round(parameters.taperlength/parameters.interval),length(inputSignal)-1);
        conv_fft = convolveWithFft([inputSignal(taperlength:-1:1,:);inputSignal;inputSignal(end:-1:end-taperlength,:)], wavelet);
        conv_fft =conv_fft(:,taperlength+1:taperlength+size(inputSignal,1));
    else
        conv_fft = convolveWithFft(inputSignal, wavelet);
    end

      
    % Create averaged power if we can
    if canAverage
        
      % Merge multiple trial
      if hasMultiTrial
        reshaped = reshape(conv_fft,r,c);
        conv_fft = mean(reshaped,2); 
      end
      
      % Calculate energy
      energy = abs(conv_fft).^2;
      
      % Save output values
      PowerMatrix(i,:) = energy;
      CoeffMatrix(i,:) = conv_fft;

      % Calc and save per trial power
      if needPerTrial
        for tr = 1 : c
          thisTrialPower       = PowerMatrixTrial{tr};
          thisTrialPower(i,:)  = abs(reshaped(:, tr)').^2;
          PowerMatrixTrial{tr} = thisTrialPower;
        end
      end
         
    else
        
      % Collect individual coefficients
      if hasMultiTrial
        reshaped = reshape(conv_fft,r,c)';
        conv_fft = reshaped';
      
         
        % Calc and save per trial power
        for tr = 1 : c
          currentMatrix        = PowerMatrixTrial{tr}; 
          currentMatrix(i,:)   = abs(conv_fft(tr,:)).^2;
          if hasMultiTrial
            PowerMatrixTrial{tr} = currentMatrix;
          end
        end
      else
        PowerMatrix(i,:) = abs(conv_fft).^2;
        CoeffMatrix(i,:) = conv_fft;
      end
      
    end
  end
  %% -------------------------
  
%   multiWaitbar(lastTitle, 'Close');
  
end