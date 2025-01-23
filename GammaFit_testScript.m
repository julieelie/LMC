% This code generates some inter spike arrival times that follow a time
% varying rate and a gamma noise as defined by the parameters of the gamma
% function. The rate at each spike is defined by the scale and shape of the
% gamma as 1./(ScaleData.*Shape).
% case study of an input signal with
FS = 1000; % sampling frequency of 1000Hz
TR = 10; %Time resolution in ms set to 10ms (work on a case where there is coherence up to 1000/(TR*2) Hz ->50Hz

% generating N random numbers from a normal distribution of mean Mu and
% standard deviation Sigma
N=1000;
Mu = 5;
Sigma = 2;
rng('shuffle')
XSamp_raw = normrnd(Mu, Sigma, N,1);

% filter the regressor so that it gets a time structure (the
% corresponding frequency of a time resolution TR where TR is 2*Tau in the
% convolution or 2*sd, F = 1/(2*Tau*pi) = 1/(TR*pi) with TR in Hz
[z,p,k] = butter(6,(1000/(TR*pi))/(FS/2),'low');
sos_band = zp2sos(z,p,k);
XSamp = filtfilt(sos_band,1,XSamp_raw);

figure(1);clf;  hold on; plot(XSamp_raw, 'k-', 'LineWidth',2);plot(XSamp, 'r-', 'LineWidth',2); legend('XSAmp_raw', 'XSamp')

% let's generate some scale data with the link function being the log instead of
% the inverse function which is the canonical function for a Gamma
% distribution. Pb with inverse fucntion is that it is illdefined for some
% particular values
% let's generate scale data with the gamma function for input of our Gamma
% GLM

B0 = -5;
B1 = 1; %0.2
Shape = 8; % The shape is fixed for all samples, only scale varies
ScaleData = exp(B0 + B1.*XSamp);
% Now add the gamma noise to the scaled data
Ydata = gamrnd(Shape, ScaleData);

figure(2); clf; histogram(Ydata, 'BinWidth',1); xlabel('Ydata: inter event intervals')

% Go from ISI (Ydata) to spike patterns
TimeFirstSpike = 5;
Duration = round(sum(Ydata) + TimeFirstSpike);
SpikePattern = cell(1,1);
SpikePattern{1} = zeros(1,Duration);
SAT = round([Ydata(1); cumsum(Ydata)]+TimeFirstSpike);
SpikePattern{1}(SAT) = 1;
TRFactors = 1:10;
Gamma  = nan(length(TRFactors),2);
SSE_Gaussconv = nan(length(TRFactors),1);
SSE_SSKernel = nan(length(TRFactors),1);
for tt = 1:length(TRFactors) %TR = 100; %window size in ms
    % Now get the time varying rate
    % Gaussian window of 2*std equal to TR (68% of Gaussian centered in TR)
    nStd =Duration/10; % before set as 4
    Tau = TRFactors(tt)*TR/2;
    T_pts = (0:2*nStd*Tau) - nStd*Tau; % centered tpoints around the mean = 0 and take data into account up to nstd away on each side
    Expwav = exp(-0.5*(T_pts).^2./Tau^2)/(Tau*(2*pi)^0.5);
    Expwav = Expwav./sum(Expwav);
    % convolution
    YRate = cell(1,1);
    YRate{1} = conv(SpikePattern{1}, Expwav,'same');
    % Instead of a convolution with a Gaussian use a kernel density estimate of the rate 
    t = 1:Duration;
    [KDE,t,Error]=kde_wrapper(SAT,t,1);
    % Now try a time varying kernel for rate estimation (nearest neighbor?)
    NNKE = zeros(length(SAT),Duration);
    ONNK = TRFactors(tt);% Order of the nearest neighbor kernel
    for ss=1:length(SAT)
        NNKE(ss,SAT(ss))=1;
        if (ss-ONNK)>0 && (ss+ONNK)<length(SAT)
            Tau = abs(SAT(ss)-min(SAT(ss+ONNK), SAT(ss-ONNK)));
        elseif (ss-ONNK)>0
            Tau = abs(SAT(ss)-SAT(ss-ONNK));
        elseif (ss+ONNK)<length(SAT)
            Tau = abs(SAT(ss)-SAT(ss+ONNK));
        else
            warning('ISSUE CANNOT FIND THE %dth nearest spike!!', ONNK);
            keyboard
        end
        Tau = max(Tau,TR/2);
        % construct the Gaussian
        T_pts = (0:2*nStd*Tau) - nStd*Tau; % centered tpoints around the mean = 0 and take data into account up to nstd away on each side
        Expwav = exp(-0.5*(T_pts).^2./Tau^2)/(Tau*(2*pi)^0.5);
        Expwav = Expwav./sum(Expwav);
        NNKE(ss,:) = conv(NNKE(ss,:), Expwav,'same');
    end
    NNKE = sum(NNKE);
    
    % This is the expected rate estimated at every interspike interval
    % (ExpectedRatet = SAT;)
    ExpectedRate = 1./(ScaleData.*Shape);
    LimTime = 1000; % for plotting the spike pattern just take the first 1000 time point (first second)
    figure(3); clf; subplot(4,1,1); hist(SAT(SAT<=LimTime),LimTime); hold on; yyaxis right; plot(SAT(SAT<=LimTime),ExpectedRate(SAT<=LimTime), '-k', 'LineWidth',2)
    subplot(4,1,2); hist(find(SpikePattern{1}(1:LimTime)),LimTime); hold on; yyaxis right; plot(YRate{1}(1:LimTime), '-r', 'LineWidth',2)
    % Correct the time of spikes by the time varying rate to see what is
    % the underlying homogenous Gamma/Poisson process, we should at least
    % retrieve the shape
    [ISIPerStim_Rescaled, ISIvecold, YPatterns_Rescaled]=time_rescaled_ISI(SpikePattern, NNKE);
    figure(3);subplot(4,1,3); hist(find(YPatterns_Rescaled{1}(1:LimTime)),LimTime);
    ExpRate = resample(ExpectedRate, SAT(2:end),1);
    figure(3); subplot(4,1,4);plot(SAT(SAT<=LimTime),ExpectedRate(SAT<=LimTime), '-k', 'LineWidth',2);hold on;plot(ExpRate(1:LimTime), '--k', 'LineWidth',2);hold on; yyaxis left; plot(YRate{1}(1:LimTime), '-r', 'LineWidth',2); hold on; plot(NNKE(1:LimTime), '-b', 'LineWidth',2); legend({'Expected Rate' 'Rescaled Expected Rate' 'calculated Rate Gaussian conv' 'calculated NN kernel'}); ylim([0 0.4])
    
    % suplabel(sprintf('SSE Gaussian convolution = %.2f SSE NNKE = %.2f', sum((YRate{1}(1:length(ExpRate))-ExpRate').^2), sum((NNKE(1:length(ExpRate))-ExpRate').^2)),'t')
    suplabel(sprintf('SSE Gaussian convolution = %.2f SSE NNKE = %.2f', sum((YRate{1}(1:length(NNKE))-ExpRate(1:length(NNKE))').^2), sum((NNKE(1:length(NNKE))-ExpRate(1:length(NNKE))').^2)),'t')
    try 
        SSE_Gaussconv(tt) = sum((YRate{1}(1:length(ExpRate))-ExpRate').^2);
    catch
        SSE_Gaussconv(tt) = sum((YRate{1}(1:length(YRate{1}))-ExpRate(1:length(YRate{1}))').^2);
    end
    try
        SSE_SSKernel(tt) = sum((NNKE(1:length(ExpRate))-ExpRate').^2);
    catch
        SSE_SSKernel(tt) = sum((NNKE(1:length(YRate{1}))-ExpRate(1:length(YRate{1}))').^2);
    end
    
    FIG = figure(6)
    clf
    subplot(4,1,1)
    scatter(XSamp,ScaleData,  'filled'); ylabel('Scale Data '); xlabel('XSamp')
    subplot(4,1,2)
    scatter(XSamp, Ydata, 'filled'); xlabel('XAmp'); ylabel('Gamma Y Data')
    hold on
    scatter(XSamp, ScaleData, 'filled', 'MarkerFaceColor','r', 'MarkerEdgeColor','r'); 
    hold on
    scatter(XSamp, ScaleData.*Shape, 'filled', 'MarkerFaceColor','k', 'MarkerEdgeColor','k'); legend('Scaled Data + Gamma noise', 'Scaled Data', 'Original Mean')
    
    subplot(4,1,3)
    H=histogram(Ydata, 'BinWidth',1); xlabel('Ydata (ms)');
    % fit the distribution of Ydata with a gamma distribution
    Gam = fitdist(Ydata, 'Gamma'); % This fit assumes a constant shape and scale, not a changing scale
    % fit the distribution of Ydata with an exponential (Poisson process)
    Poi = fitdist(Ydata, 'Exponential');
    % add the 2 fit to the figure
    XISI = H.BinEdges(1):(H.BinEdges(end)-1);
    Xexp = pdf('Exponential',XISI,Poi.mu);
    Xgam = gampdf(XISI, Gam.a, Gam.b);
    subplot(4,1,3);cla
    yyaxis right; plot(XISI, Xexp, '-', 'LineWidth',2, 'Color', [0.929 0.694 0.125]);
    hold on; yyaxis right; plot(XISI, Xgam, '-', 'LineWidth',2);
    hold on; yyaxis left; H=histogram(Ydata, 'BinWidth',1); xlabel('Ydata (ms)');hold off
    legend({'Ydata','Exp fit', 'Gamma fit'})
    title(sprintf('Ydata fit Gamma shape parameter = %.2f Gamma scale parameter = %.2f ', Gam.a,Gam.b))
    
    subplot(4,1,4)
    H=histogram(ISIPerStim_Rescaled{1}, 'BinWidth',1); xlabel('Ydata (ms)');
    % fit the distribution of Ydata with a gamma distribution
    Gam = fitdist(ISIPerStim_Rescaled{1}', 'Gamma'); % This fit assumes a constant shape and scale, not a changing scale
    Gamma(tt,1) = Gam.a;
    Gamma(tt,2) = Gam.b;
    
    % fit the distribution of Ydata with an exponential (Poisson process)
    Poi = fitdist(ISIPerStim_Rescaled{1}', 'Exponential');
    % add the 2 fit to the figure
    XISI = H.BinEdges(1):(H.BinEdges(end)-1);
    Xexp = pdf('Exponential',XISI,Poi.mu);
    Xgam = gampdf(XISI, Gam.a, Gam.b);
    subplot(4,1,4);cla
    yyaxis right; plot(XISI, Xexp, '-', 'LineWidth',2, 'Color', [0.929 0.694 0.125]);
    hold on; yyaxis right; plot(XISI, Xgam, '-', 'LineWidth',2);
    hold on; yyaxis left; H=histogram( ISIPerStim_Rescaled{1}, 'BinWidth',1); xlabel('Ydata rescaled (ms)');hold off
    legend({'Ydata','Exp fit', 'Gamma fit'})
    title(sprintf('Ydata fit Gamma shape parameter = %.2f Gamma scale parameter = %.2f ', Gam.a,Gam.b))
    pause()
end

% Here we see that TR has a strong efect, if the estimation of the rate is
% not correct and influence too much by local noise, then estimation of
% gamma parameters using the scaling method is wrong.
figure()
subplot(2,1,1)
yyaxis left
plot(TRFactors,Gamma(:,1),'b-', 'LineWidth',2)
ylabel('Shape parameter')
xlabel('TRFactors')
yyaxis right
plot(TRFactors,Gamma(:,2),'r-', 'LineWidth',2)
ylabel('Scale parameter')
subplot(2,1,2)
plot(TRFactors, SSE_Gaussconv); hold on; plot(TRFactors, SSE_SSKernel);ylabel('SSE'); xlabel('TR Factors'); legend({'Gaussian convolution' 'SSkernel'})
% THE 4th NEAREST NEIGHBOURG KERNEL ESTIAMTE OF THE RATE SEEMS TO GIVE THE
% BEST VALUES 

% Now let's fit the data with a GLM gamma
[B_GLM,Dev, Stats] = glmfit(XSamp,Ydata,'gamma', 'link', 'log');
GamModel = fitglm(XSamp,Ydata,'linear','Distribution','gamma', 'DispersionFlag', true,'Link', 'log');

 % estimated shape
fprintf(1,'The model estimates a shape of %.1f when it was set as %.1f\n', 1/GamModel.Dispersion, Shape)
fprintf(1, 'The model estimates the intercept to be %.1f when it was set as %.1f\n', GamModel.Coefficients.Estimate(1) - log(1/GamModel.Dispersion), B0)
fprintf(1, 'The model estimates the slope to be %.1f when it was set as %.1f\n', GamModel.Coefficients.Estimate(2), B1)

