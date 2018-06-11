% createMorletWavelet creates a Morlet wavelet with base frequency
% given in parameter.
%
% wavelet = createMorletWavelet(baseFrequency, parameters)
%   f - base frequency for the wavelet
%   parameters
%           interval - sampling interval for the wavelet
%           wavenumber - if it is given, this will be the wavelet number,
%           else the default 7 will be used
function wavelet = createMorletWavelet(baseFrequency, parameters)

    %If we want to visualize the resulting wavelet, set this to 1
    PLOT = 0;
    
    % Default wave number will be 7
    if ~isfield(parameters,'wavenumber')
        wavenumber = 7;
    else
        wavenumber = parameters.wavenumber;
    end
    
    interval = parameters.interval;

    sigma_f = round(baseFrequency/wavenumber);
    sigma_t = 1/(2*pi*sigma_f);
    sigma_tsq = sigma_t^2;
    A = sqrt(1/(sigma_t * sqrt(pi)));
        
    t = -1:interval:1;
    
    if mod(length(t),2)==0
       t = [ t , 1 ];  
    end
    
    t_sq = t.^2;
    wavelet = A .* exp(-t_sq ./ (2 * sigma_tsq)) .* exp(2*1i*pi*baseFrequency*t);
    
    if PLOT==1
        figure; plot(t, real(wavelet)); title('Wavelet');
    end
    
end