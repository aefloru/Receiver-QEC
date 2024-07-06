function [corrected_signal] = correctIQImbalance()

    [I_signal, Q_signal] = RandomWaveforms();
    % Apply low-pass filtering
    fbw = 0.2; % bandwidth of the LPF

    N = length(I_signal); % number of samples
    % sample indices, from transmitter
    n = (0:(N-1))';

    % Step 2: Compute βI and βQ (DC offsets)
    beta_I = mean(I_signal);
    beta_Q = mean(Q_signal);

    % Step 3: Remove the DC offsets
    I_error = I_signal - beta_I;
    Q_error = Q_signal - beta_Q;

    %cross_correlationIQ = (alpha*sin(psi) * I_signal.^2);
   
    % Step 4: Compute α (amplitude error) g_q 
    alpha = sqrt(mean(I_error.^2) / mean(Q_error.^2));

    % Step 5: Compute sin(ψ) (phase error) phi_q
    psi = asin((mean(I_error.*Q_error))/ sqrt(mean(I_error.^2) .* mean(Q_error.^2)));

    % Step 7: Compute the correction matrix parameters
    A = 1 / alpha;
    C = -sin(psi) / (alpha * cos(psi));
    D = 1 / cos(psi);
    % Step 8: Apply the correction
    corrected_signal = zeros(2, N);
    corrected_signal(1, :) = A * (I_error);  % Corrected I, first row, second term goes to 0
    corrected_signal(2, :) = C * (I_error) + D * (Q_error);  % Corrected Q, second row

    I_corr = corrected_signal(1, :);
    Q_corr = corrected_signal(2, :); 
    y_corrected = I_corr + 1i*Q_corr; %returning complex signal 
    
    figure
%    subplot(2,1,1);
%    PlotPsd(I_signal + 1i*Q_signal);
%    title('Uncorrected Signal');

%    subplot(2,1,2)
%    PlotPsd(y_corrected);
%    title('Corrected Signal');  

figure;
% Plot the power spectral density (PSD) of the uncorrected signal in red
PlotPsd(I_signal + 1i*Q_signal, 1, "red");
hold on; % Hold the current plot to overlay the next plot
% Plot the PSD of the corrected signal in blue


PlotPsd(y_corrected, 1, "b");

% Add title and legend
title('Uncorrected and Corrected Signals');
legend('Uncorrected Signal with Images', 'Corrected Signal');
hold off; % Release the plot hold


end

 function [real_val, img_val] = DirectDownConversionDemo()
    % Define parameters
    N = 1e5; % number of samples
    n = (0:N-1)'; % sample indices
    f_c = 0.1; % normalized carrier frequency
    phi_c = 0; % carrier phase
    phi_q = 0.0; % quadrature phase error
    g_q = 0.3; % quadrature gain imbalance

    % Generate a complex bandpass signal
    f_x = 0.05; % signal frequency
    x = exp(1i*2*pi*f_x*n); % complex exponential at f_x
    % 
    % % Apply quadrature error
    g1 = (1 + g_q * exp(-1i*phi_q)) / 2;
    g2 = (1 - g_q * exp(1i*phi_q)) / 2;
    y_bb = x * exp(-1i*phi_c) .* g1 + conj(x) * exp(-1i*phi_c) .* g2; % Baseband signal with quadrature error
    
    % Apply low-pass filtering
    fbw = 0.2; % bandwidth of the LPF

    % Sample the filtered signal (ADC)
    y_adc = y_bb; % 
    real_val = real(y_adc);
    img_val = imag(y_adc);
    % Plot the spectrum of the baseband signal
    plot_ft(y_adc, 'Baseband Signal Spectrum with Quadrature Error');

end

function y = ApplyLpf(x, fbw, doPlot)
    
    % build filter
    M = 100;
    n = ((-M):(M))';
    w0 = 2*pi*fbw/2;
    h = sin(w0*n)./(pi*n);
    h(n == 0) = w0/pi;
    h = h .* kaiser(length(h), 10);
    
    % apply filter
    y = conv(x, h, 'valid');
    
    % We are simulating an analog filter, so there will still be
    % noise sources before we sample at the ADC
    % Let's add some noise back in to make the plots look nice.
    % (otherwise, the noise floor would get shaped by the filter too)
    y = y + 10^(-40/20)*randn(length(y), 1);
    
    % plot
    if ( (nargin >= 3) && doPlot )
        N = length(x);
        f = (-1/2):(1/N):(1/2-1/N);        
        wx = blackman(length(x));
        wy = blackman(length(y));
        X = fftshift(fft(wx.*x, N)) / sum(wx);
        Y = fftshift(fft(wy.*y, N)) / sum(wy);
        H = fftshift(fft(h, N));        
        figure;        
        subplot(2,1,1);
        hold on;
        plot(f, 20*log10(abs(X)), 'b');
        plot(f, 20*log10(abs(H)), 'r');
        legend('X(f)', 'H(f)');
        title('Input X(f)');
        xlabel('Frequency (normalized)');
        ylabel('Amplitude (dB)');
        xlim([-0.5 0.5]);
        ylim([-150 0]);
        subplot(2,1,2);
        hold on;
        plot(f, 20*log10(abs(Y)), 'b');
        plot(f, 20*log10(abs(H)), 'r');
        legend('Y(f)', 'H(f)');
        title('Output Y(f)');
        xlabel('Frequency (normalized)');
        ylabel('Amplitude (dB)');
        xlim([-0.5 0.5]);
        ylim([-150 0]);
    end
    
end

function plot_ft(x, titleStr)
    N = length(x);
    f = (-N/2):(N-1)/2;
    f = f / N; % Normalize frequencies
    
    X = fftshift(fft(x));
    figure;
    plot(f, 20*log10(abs(X)/max(abs(X))));
    ylabel('Amplitude (dB)');
    xlabel('Frequency (normalized)');
    title(titleStr);
end


function PlotPsd(x, fs, color)

   % if (nargin < 2)
   %     fs = 1;
   % end
    
    [pxx, f] = GetPsd(x, fs);
    
    plot(f/1e6, 10*log10(pxx), color);
    xlabel('Frequency (MHz)');
    ylabel('PSD (dBFS/Hz)');
    
end

function [pxx, f] = GetPsd(x, fs)

    if (nargin < 2)
        fs = 1;
    end

    N = 2^floor(log2(length(x)/16));
    M = N/4;
    L = N*4;
    w = blackman(N);

    [pxx, f] = pwelch(x, w, M, L, fs, 'centered');
    
end
