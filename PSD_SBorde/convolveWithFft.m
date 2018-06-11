%   convolveWithFft performs the convolution with FFT.
%
%   CONV_FFT = convolveWithFft(SIGNAL, RESPONSE) convolves the signal
%   with the given response function. Pad both vector to the same, 2-power
%   length vector, calculates its FFT, and performs IFFT after
%   multiplication in frequency domain.
function conv_fft = convolveWithFft(signal, response)

    signal_length = length(signal);
    response_length = length(response);
    new_length = signal_length + response_length - 1;
    new_length_pow = pow2(nextpow2(new_length));
    half_response_length = round((response_length-1)/2);
        
    [r,c] = size(signal);
    if r>c
        signal=signal';
    end
    [r,c] = size(response);
    if r>c
        response=response';
    end
    
    conv_fft = ifft( fft(response, new_length_pow) .* fft(signal, new_length_pow) );
    conv_fft = conv_fft(1:new_length);
    conv_fft = conv_fft(half_response_length+1:end-half_response_length);
    
end