function [CoherencyT, fo, Coherence, Coherence_up, Coherence_low, stP, H] = multitapercoherence_JN_fast(x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers)
% used to be function [y, fo, meanP, Pupper, Plower, stP] =
% df_mtchd_JN(varargin) from STRFLab Theunissen Lab
%function [yo, fo, yJ, yupJ, ylowerJ, stJ]=mtcsg(x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers);
% Multitaper Cross-Spectral Density, jacknived estimates and errors
%only meanP(:,1,2), Pupper(:,1,2), Plower(:,1,2) is the correct jack-knifed
%coherence.
%These values are the absolute values of the coherency, to get coherence, these values
%must be squared.
%y is the original cross spectrum without jack-knifing or normalizing to get coherency.
% function A=mtcsd(x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers)
% x : input time series, a 2 column matrix where each column contains the
% variable between which you want to calculate the coherence
% nFFT = number of points of FFT to calculate (default 1024) needs to be at
% least the value of the Nyquist limit for your signals (Fs/2)
% Fs = sampling frequency (default 2)
% WinLength = length of moving window (default is nFFT)
% nOverlap = overlap between successive windows (default is WinLength/2)
% NW = time bandwidth parameter (e.g. 3 or 4), default 3
% nTapers = number of data tapers kept, default 2*NW -1
%I've changed this program to output the magnitude of coherency, or coherence.
% output yo is yo(f)
%
% If x is a multicolumn matrix, each column will be treated as a time
% series and you'll get a matrix of cross-spectra out yo(f, Ch1, Ch2)
% NB they are cross-spectra not coherences. If you want coherences use
% mtcohere
%
% Original code by Partha Mitra - modified by Ken Harris
% Also containing elements from specgram.m
% 
% Modified by Anne, 2002.
% Modified by Julie E Elie, 2021

%% Sorting inputs and setting default values
if nargin<2
    nFFT = 1024;
end

if nargin<3
    Fs = 2;
end
if nargin<4
    WinLength = nFFT;
end
if nargin<5
    nOverlap = WinLength/2;
end
if nargin<6
    NW = 3;
end
if nargin<7
    Detrend = ''; 
end
if nargin<8
    nTapers = 2*NW -1;
end


WinLength = round(WinLength);
winstep = WinLength - nOverlap;
winstep = round(winstep);
nFFT = round(nFFT);

nChannels = size(x, 2);   % In most cases nChannels will be 2 : the
                          % number of time series.
nSamples = size(x,1);

% check for column vector input
if nSamples == 1 
	x = x';
	nSamples = size(x,1);
	nChannels = 1;
end

% calculate number of FFTChunks per channel
nFFTChunks = round(((nSamples-WinLength)/winstep));
% turn this into time, using the sample frequency
% t = winstep*(0:(nFFTChunks-1))'/Fs;

% calculate Slepian sequences.  Tapers is a matrix of size [WinLength, nTapers]

%[JN,y,stP] = make_slepian(x,WinLength,NW,nTapers,nChannels,nFFTChunks,nFFT,Detrend,winstep);
% allocate memory now to avoid nasty surprises later
% stP = zeros(nFFT, nChannels, nChannels);
varP = zeros(nChannels, nChannels,nFFT);
[Tapers, ~]=dpss(WinLength,NW,nTapers, 'calc');
Periodogram = complex(zeros(nFFT, nTapers, nChannels)); % intermediate FFTs
Temp1 = complex(zeros(nFFT, nTapers)); %Temps are particular psd or csd values for a frequency and taper
Temp2 = complex(zeros(nFFT, nTapers));
Temp3 = complex(zeros(nFFT, nTapers));
eJ = complex(zeros(nFFT,1));
JN = complex(zeros(nFFTChunks,nChannels, nChannels,nFFT));
%jackknifed cross-spectral-densities or csd. Note: JN(.,1,1,.) is
%the power-spectral-density of time series 1 and JN(.,2,2,.) is the
%psd of time series 2.  Half-way through this code JN(.,1,2,.)
%ceases to be the csd of 1 and 2 and becomes the abs coherency of 1
%and 2. 
y=complex(zeros(nChannels, nChannels,nFFT)); % output array for csd
Py=zeros(nChannels, nChannels, nFFT); % output array for psd's


% New super duper vectorized alogirthm
% compute tapered periodogram with FFT 
% This involves lots of wrangling with multidimensional arrays.

TaperingArray = repmat(Tapers, [1 1 nChannels]);
for jj=1:nFFTChunks
	Segment = x((jj-1)*winstep + (1:WinLength), :);
    if (~isempty(Detrend))
        Segment = detrend(Segment, Detrend);
    end
	SegmentsArray = permute(repmat(Segment, [1 1 nTapers]), [1 3 2]);
	TaperedSegments = TaperingArray .* SegmentsArray;

	Periodogram(:,:,:) = fft(TaperedSegments,nFFT);

	% Now make cross-products of them to fill cross-spectrum matrix
	for Ch1 = 1:nChannels
		for Ch2 = Ch1:nChannels % don't compute cross-spectra twice
			Temp1 = squeeze(Periodogram(:,:,Ch1));
			Temp2 = squeeze(Periodogram(:,:,Ch2));	
			Temp2 = conj(Temp2);
			Temp3 = Temp1 .* Temp2;

			%eJ and eJ2 are the sum over all the tapers.
			eJ=sum(Temp3, 2)/nTapers;
			JN(jj,Ch1, Ch2,:) = eJ;  % Here it is just the
						% cross-power for one
						% particular chunk.
			y(Ch1, Ch2,:)= y(Ch1,Ch2,:) + shiftdim(eJ,-2);  % y is
                                                           % the
                                                           % sum of
                                                           % the cross-power
		end
	end
end 

% now fill other half of matrix with complex conjugate
for Ch1 = 1:nChannels
	for Ch2 = (Ch1+1):nChannels % don't compute cross-spectra twice
		y(Ch2, Ch1,:) = y(Ch1,Ch2,:);
        Py(Ch1, Ch2, :) = atanh(abs(y(Ch1,Ch2,:)./sqrt(abs(y(Ch1,Ch1,:)).*abs(y(Ch2,Ch2,:)))));
        if (Ch1 == 1) && (Ch2==2)
            H = y(Ch1, Ch2,:)./y(Ch1,Ch1,:);
        end
	end
end

% Set NaN values to zero -- they are the result of both numerator
% and denominator in the previous step being 0. The cross-spectral
% power in this case is zero.
nans = isnan(Py(Ch1, Ch2,:));
if any(nans)
    Py(Ch1, Ch2,nans) = 0;
end


for jj = 1:nFFTChunks
    JN(jj,:, :, :) = abs(y - squeeze(JN(jj,:, :,:))); % This is where
                                                    % it becomes
                                                    % the JN
                                                    % quantity
                                                    % (the delete one)
    for Ch1 = 1:nChannels
        for Ch2 = (Ch1+1):nChannels  
            % Calculate the transformed coherence
            %JN(j, Ch1, Ch2,:) =atanh(abs(JN(j,Ch1,Ch2,:))./sqrt(abs(JN(j,Ch1,Ch1,:)).*abs(JN(j,Ch2,Ch2,:))));
            JN(jj, Ch1, Ch2,:) =atanh(real(JN(jj,Ch1,Ch2,:))./sqrt(abs(JN(jj,Ch1,Ch1,:)).*abs(JN(jj,Ch2,Ch2,:))));
            % Obtain the pseudo values
            JN(jj,Ch1, Ch2,:) = shiftdim(nFFTChunks*Py(Ch1, Ch2,:),-1) - (nFFTChunks-1)*JN(jj, Ch1,Ch2,:);
        end
    end  
end

% Set NaN values to zero -- they are the result of both numerator
% and denominator in the previous step being 0. The cross-spectral
% power in this case is zero.
nans = isnan(JN);
if any(nans(:))
    JN(nans) = 0;
end

meanP=squeeze(mean(JN,1));
for Ch1=1:nChannels
    for Ch2=Ch1:nChannels
        varP(Ch1, Ch2,:) = (1/nFFTChunks)*var(squeeze(JN(:,Ch1, Ch2,:)),1);
    end
end

%upper and lower bounds will be 2 standard deviations away.
stP=sqrt(varP);

size(stP); % Should be NChannels, NChannels,nFFT
size(meanP);% Should be NChannels, NChannels,nFFT

Pupper = meanP + 2*stP;
Plower = meanP - 2*stP;
meanP = tanh(meanP);
Pupper = tanh(Pupper);
Plower = tanh(Plower);

		
% set up f array
if ~any(any(imag(x)))    % x purely real
	if rem(nFFT,2)  % nfft odd
		select = 1:(nFFT+1)/2;
	else
		select = 1:nFFT/2+1;
	end
    meanP = meanP(:,:, select);
    Pupper = Pupper(:,:,select);
    Plower = Plower(:,:,select);
% 	y = y(select,:,:);
else
	select = 1:nFFT;
end
	
fo = (select - 1)'*Fs/nFFT;

% we've now done the computation.

% Output the coherency in Time (non-jackknife) and the coherence
Coherence = squeeze(meanP(1,2,:));
Coherence_up = squeeze(Pupper(1,2,:));
Coherence_low = squeeze(Plower(1,2,:));
% y contains the cross and auto correlation spectral densities of the 2
% signals in x so the time coherency is obtained by:
CoherencyT = squeeze(ifft(y(1,2,:)./(abs(y(1,1,:)) .* abs(y(2,2,:))).^0.5));

% the rest of this code is stolen from
% specgram and just deals with the output stage

if nargout == 0
    % take abs, and plot results
    newplot;
    for Ch1=1:nChannels
        for Ch2 = 1:nChannels
            subplot(nChannels, nChannels, Ch1 + (Ch2-1)*nChannels);
            plot(f,20*log10(abs(y(Ch1,Ch2,:))+eps));
            grid on;
            if(Ch1==Ch2) 
                ylabel('psd (dB)'); 
            else 
                ylabel('csd (dB)'); 
            end
            xlabel('Frequency');
        end
    end
end
% Original mtcsd Written by Kenneth D. Harris
% Jackknife estimates over FFT chunks by Anne Hsu and Frederic Theunissen
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu
