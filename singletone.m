 function [ I_signal , Q_signal ] = DirectDownConversionDemo()
% Define parameters
 N = 1e5 ; % number of samples
 n = (0: N -1); % sample indices
 f_c = 0.1; % normalized carrier frequency
 phi_c = 0; % carrier phase
 phi_q = 0.0; % quadrature phase error
 g_q = 0.8; % quadrature gain imbalance

 % Generate a complex bandpass signal
 f_x = 0.1; % signal frequency
 x = exp (1i *2* pi* f_x * n ) ; % complex exponential at f_x

 plot_ft(x , "Ideal Receiver Output") ;

% Apply quadrature error
g1 = (1 + g_q * exp ( -1i * phi_q ) ) / 2;
g2 = (1 - g_q * exp (1i * phi_q ) ) / 2;
y_bb = x * exp ( -1i * phi_c ) .* g1 + conj ( x ) * exp ( -1i * phi_c ) .* g2 ; % Baseband signal with quadrature error


 % Sample the filtered signal (ADC)
 y_adc = y_bb ;
 I_signal = real( y_adc ) ;
 Q_signal = imag( y_adc ) ;

 % Plot the spectrum of the baseband signal
 plot_ft( y_adc , "Received Signal with QE") ;

 end


 function plot_ft(x, t)
    N = length(x);
    % Compute frequency vector for plotting
    f = (-N/2):(N-1)/2;
    f = f / N; % Normalize frequencies
    
    figure;
    plot(f, 20*log10(abs(fftshift(fft(x)))/N));
    ylabel('Amplitude (dB)');
    xlabel('Frequency (normalized)');
    title(t);

end