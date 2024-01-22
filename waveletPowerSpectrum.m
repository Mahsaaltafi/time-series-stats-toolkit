function mPSD = waveletPowerSpectrum(In, BetaTimes)
% WAVELETPOWERSPECTRUM Compute the mean Power Spectral Density (mPSD) around beta oscillation peaks.
%
%   This function calculates the wavelet-based power spectrum of a signal 'In'
%   sampled at 1250 Hz, focusing on the periods of beta oscillations. The analysis
%   window is centered around the beta oscillation peaks with a margin of +/- 240 ms
%   (300 samples at 1250 Hz).
%
%   This code has been formated for readability by ChatGPT-4.
%
%   Inputs:
%       In - Vector of the input signal sampled at 1250 Hz.
%       BetaTimes - Matrix of beta event times with each row as [start peak end] in seconds.
%
%   Outputs:
%       mPSD - A matrix with rows corresponding to each analyzed frequency band.
%              The columns represent:
%              mPSD(:, 1) - Frequencies in Hz.
%              mPSD(:, 2) - Mean of the wavelet coefficients over the analysis window.
%              mPSD(:, 3) - Standard error of the mean (SEM) of the wavelet coefficients
%                            over the analysis window.
%
%   Example:
%       signal = randn(1, 30000); % Example signal
%       betaTimes = [0.85 1 1.15; 2.35 2.5 2.65]; % Example beta event times in seconds
%       mPSD = waveletPowerSpectrum(signal, betaTimes);
%
%   Note:
%       The BetaTimes input should have 3 columns representing the start, peak, 
%       and end times of beta oscillation events in seconds.
%
%   See also WAVELET

  % Convert BetaTimes to sample indices
    BT = round(BetaTimes(:,2) * 1250);

    % Calculate wavelet transform for each beta event
    for i = 1:length(BT)
        [P(:,:,i), F] = cwt(In(BT(i)-300:BT(i)+300), 'amor', 1250);
    end

    % Filter frequencies below 80 Hz
    myf = F < 80;

    % Prepare data for the spectrogram
    P2 = real(P(myf,:,:));
    Spectrogram.data = nanmean(P2, 3);
    Spectrogram.freq = F(myf);
    Spectrogram.time = [-300:300] / 1.25;

    % Plot the spectrogram
    h = pcolor(Spectrogram.time, Spectrogram.freq, Spectrogram.data);
    set(h, 'EdgeColor', 'none');
    set(gca, 'layer', 'top');
    shading INTERP;
    box off;
    set(gca, 'TickDir', 'out');
    xlim([-200 200]);
    yticks([10 15 20 25 30 35 40 45 50 55 60 65 70 75 80]);
    colormap jet;
    clf;

    % Calculate mean power spectrum density
    mPSD(:,1) = F(myf);
    PSD = abs(squeeze(real(P(myf,301,:))));
    mPSD(:,2) = mean(PSD, 2);
    mPSD(:,3) = sem(PSD');
    mPSD = sortrows(mPSD, 1);
    clf;

    % Plot the power spectrum
    errorbar(mPSD(:,1), mPSD(:,2), mPSD(:,3));
end


function Plot_powerspectrum(mSpectrum) 
    % Plot function for power spectrum
    figure;
    boundedline(mSpectrum(:,1), mSpectrum(:,2), mSpectrum(:,3));
    xticks(0:5:80);
    view(90, -90);
end

function SEMean=sem(A)
    % Computes standard error of mean
    if sqrt(length(A(find(isnan(A)==0))))==0
       SEMean=NaN;
    else
       V=isnan(A)==0;
       SEMean=nanstd(A)./sqrt(sum(V,1));
    end
end

