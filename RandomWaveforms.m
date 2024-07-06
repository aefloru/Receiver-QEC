function [I, Q] = RandomWaveforms()

    % The previous MatlabDemo generated sinusoidal signals and discussed
    % frequency-domain analysis via FFTs.  In this new demo, we consider
    % signals that are more closely related to modern communications
    % signals.  We could work with real communications signals, but we
    % might find it easier for our purposes to approximate communications
    % signals with Gaussian noise.  When we deal with "random" signals
    % like Gaussian noise, then straight FFTs will not give us the best
    % estimate of the signal's Power Spectral Density (PSD).  Instead,
    % there are other techniques for estimating PSD.
    %
    % Just like with the FFT, for this course you do not need to fully
    % understand what is happening under-the-hood.  We will simply use
    % our PSD estimate as a visualization tool.
    
    % number of samples to generate
    N = 1e4;
    
    % sampling rate in Hz
    FS = 100e6;
    
    % sample indices
    n = (0:(N-1))';
    
    % generate a narrowband modulated waveform
    F = [20e6 -30e6 0e6];       % frequency
    P = 10.^([-12 -20 -25]/10); % power
    B = [0e6 5e6 10e6];         % bandwidth (approximate)
    x = zeros(N,1);
    for i = 1:length(P)
        if (B(i) == 0)
            % special case: generate CW
            x0 = sqrt(P(i))*exp(1i*2*pi*rand(1,1))*ones(N,1);
        else
            % random narrowband waveform
            U = ceil(FS/B(i));
            M = ceil(N/U);
            x0 = sqrt(P(i))*sqrt(2)/2*(randn(M,1) + 1i*randn(M,1));
            x0 = resample(x0, U, 1);
        end
        x0 = x0(1:N) .* exp(1i*2*pi*F(i)/FS*n);
        x = x + x0;
    end    
    
    % add some noise just to make the plots look more realistic    
    x = x + 10^(-50/20)*sqrt(2)/2*(randn(N,1) + 1i*randn(N,1));       
    
    % apply baseband-equivalent receiver model to create image
    gq =  0.98;
    pq = -0.02;
    g1 = (1/2)*(1 + gq*cos(pq) - 1i*gq*sin(pq));
    g2 = (1/2)*(1 - gq*cos(pq) - 1i*gq*sin(pq));
    y = g1*x + g2*conj(x);
    I = real(y);
    Q = imag(y);
    % plot the evaluate the signals in the frequency domain
    % use the PlotPsd() function provided below    
    
    figure;
    %set(gcf, 'WindowStyle', 'docked');
    subplot(2,1,1);
    PlotPsd(x, FS);    
    title('Original Signal');%title('|X(f)|^2');
    %subplot(2,1,2);
    %PlotPsd(y, FS);    
    %title('|Y(f)|^2');

end

function PlotPsd(x, fs)

    if (nargin < 2)
        fs = 1;
    end
    
    [pxx, f] = GetPsd(x, fs);
    
    plot(f/1e6, 10*log10(pxx));
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