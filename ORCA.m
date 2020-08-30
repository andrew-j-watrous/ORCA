function [NormAmp,bandfreqs,bands,bandamps,bandphases,recon_stats,params] = ORCA(signal,params)
%   Oscillatory Reconstruction Algorithm (ORCA)
%   Andrew J Watrous, Ph.D.
%   Preprint posted on BiorXiv and paper submitted for publication:
%   https://github.com/andrew-j-watrous/ORCA

%   Dependencies: 
%   1) Wavelet Toolbox in Matlab (for cwt.m & icwt.m)
%   2) "EEG Toolbox" from Kahana lab.
%   See: http://memory.psych.upenn.edu/Software                
%   3) 'linspecer' from Matlab FEX if you wish to plot bands

%   Inputs
%   Requires signal and params
%
%   Signal is the signal you wish to analyze. 
%
%
%   Params must include the following two fields:
%   params.freqrange:        [min max] Broadband frequency range to consider (Hz)
%
%   params.srate:       Original sampling rate of the input signal(Hz) (e.g. 1000 Hz
%                       sampling: params.srate = 1000). If you elect to
%                       downsample the signal using params.downsample=1,
%                       you should *not* change this number to your desired
%                       sampling rate but instead always define this as the
%                       original sampling rate of the input signal.
%
%   params.bandMethods:       Cell array listing which methods to try, options:                    
%                              Example: params.bandMethods = {[1 4;4 8;8 12;12 30;30 100],'SCV','CoD','PeakPick'};
%
%                      'CoD'.     Define bands based on local minima in
%                                 coefficient of determination (r2) across
%                                 frequencies
%                                 
%                      'SCV'      Define bands based on local minima in
%                                 spectral Coefficient of variation across
%                                 frequencies
%
%                      'PeakPick' Old 1/f fitting method used in MODAL
%                      (Watrous et al., 2018 eLife). This method places a
%                      band boundary at frequencies transitioning above/below
%                      the 1/f "background" spectrum.
%                      
%                      [N X 2]    User defined band edges using N bands. 
%                                 Min and max freq are set by .5 and Nyquist frequency
%                      
%


%   In addition to the required inputs above, there are many optional input parameters:
%   Most of these are Booleans/logical

%   params.downsample        Boolean. Downsample signal (if possible, based
%                            on maximum frequency of interest). This will
%                            speed computation, sometimes quite
%                            considerably.  Default: 0/no. 
%                            
%   params.broadband_filter: Boolean. Wideband filter signal to freqrange (default: 0/no)
%
%   params.min_bw:           minimum width of a band (Hz; default .5).  Too narrow of
%                            bands leads to filtering errors, particularly at very low frequencies
%
%   params.verbose           Boolean for verbose mode (default: on/1)
%
%   params.crossval_trainIDX Boolean vector of length equal to signal that defines which signal samples 
%                            are to be used for training during crossvalidation. Remaining signal samples                            
%                            are used for testing. By default, the first
%                            half of the signal is used if not-specified by
%                            user. Note that particular windows, such as
%                            pre-trial baselines, can be used to define
%                            bands.

%   params.numControlPerms   Default: 0: Number of times to calculate r^2 using shuffled
%                            (circshifted over time axis) amplitude, phase, and frequency estimates.
%
%   params.bad_data          Boolean vector of bad data. 1 == bad data to exclude from
%                            calculations.  Must be length of signal.
                             %currently only excludes data when using 'SCV'
                             %or 'PeakPick' methods
                             
%   params.plot_bands        Boolean. Make plots of bands.

%   params.normalize_amp     Default: 1. Normalize halfwave amplitudes by
%                            their detected instantaneous frequency. This 
%                            ensures that slower activity in a band is not
%                            determined to be larger as a mere
%                            consequence of the 1/f nature of EEG.

%   params.saveFiltered      default: 0/no: save out the filtered traces
%                            from the best method in recon_stats
    


%   Outputs
%   NormAmp:                 Bands X Samples. Normalized Amplitude
%   bandfreqs                Bands X Samples. Frequency estimates using modified "frequency sliding" from MX Cohen J Neuro 2014
%   bands:                   band edges [low_freq high_freq] per band
%   bandamps :               Band X samples amplitudes
%   bandphases:              Band X samples phases (radians).
%   recon_stats:             structure containing information on
%                            reconstruction, including reconstructed signal, best_r2, best method,
%                            and band importance
%   params                   params updated with what was computed
   
    
%display header
disp('************************')


%Establish parameters/inputs
if ~isfield(params,'min_bw') %mininmum bandwidth
    params.min_bw = .1;
end

if ~isfield(params,'broadband_filter') %broadband filter, default is to broadband filter
    params.broadband_filter = 1;
end

if ~isfield(params,'saveFiltered') %broadband filter, default is to broadband filter
    params.saveFiltered = 0;
end


if ~isfield(params,'downsample')
    params.downsample = 0; %default is to not downsample
end

if ~isfield(params,'verbose')
    params.verbose = 1; %default is verbose on. Will print bands and r^2
end

%establish frequencies of interest (default: from .5 to nyquist/1.15 based on Frequency sliding algo)
nyq = (params.srate/2)-.0001; %make it slightly less than nyquist to ensure
maxfreq = floor((nyq/1.15)); %based on frequency sliding trans width of .15
minfreq = params.freqrange(1);
if isfield(params,'freqrange')
    if params.freqrange(2)>maxfreq %if user inputs maxfreq which is too high based on frequency sliding algo
        params.freqrange(2) = maxfreq;
        disp(['Updating Maximum Frequency to ',num2str(maxfreq)]);
    end
    
    params.wavefreqs = logspace(log10(params.freqrange(1)),log10(params.freqrange(2)),200); %logspace from min to max freq
else
    %we don't screw up broadband filtering
    expp1 = log10(minfreq);expp2 = log10(maxfreq);
    params.wavefreqs = logspace(expp1,expp2,200); %logspace 1-100 Hz
end

%crossvalidation
if ~isfield(params,'crossval_trainIDX')
    params.crossval_trainIDX = zeros(1,length(signal));
    params.crossval_trainIDX(1:round(length(signal)/2))=1; %default is to use first half of signal to generate bands
end

%by default, don't calculate r2 using circshifted amp, phase,freq
if ~isfield(params,'numControlPerms')
    params.numControlPerms=0;
elseif isfield(params,'numControlPerms') && params.numControlPerms>0
    disp(['Will calculate r^2 using ',num2str(params.numControlPerms),' shuffled estimates'])
end

if ~isfield(params,'plot_bands')
    params.plot_bands=0; %don't plot by default
end


if ~isfield(params,'normalize_amp') || params.normalize_amp
    params.normalize_amp=1;
    disp('Normalizing wave amplitude by frequency')
elseif isfield(params,'normalize_amp') && ~params.normalize_amp
    disp('Skipping wave amplitude normalization by frequency')
end


%%%%%%%%% Above is defining parameters
%%%%%%%%% Below is computational steps

%make sure signal is the correct shape
if size(signal,2)>1
    signal = signal';
end


%downsample signal to save computational time if our max frequency of interest is < 1/3
%of original sampling rate
params.downsample_factor = floor((params.srate/max(params.wavefreqs))/3); %how much to downsample signal
if params.downsample_factor>1 && params.downsample
    signal = downsample(signal,params.downsample_factor);
    params.srate = params.srate/params.downsample_factor;
    disp(['Downsampling signal to ',num2str(params.srate),' Hz to save time'])
else
    disp('Downsampling is OFF')
end


%ensure the signal is mean centered so that the
%hilbert transform and pow/phase estimates are valid
signal= signal-nanmean(signal);

%broadband filter the signal to within the broadest frequency range we
%might care about
if params.broadband_filter
    
    disp('Broadband filtering is ON')
    signal = buttfilt(signal,[min(params.wavefreqs) max(params.wavefreqs)],...
        params.srate,'bandpass',4);
else
    disp('Broadband filtering is OFF')
    
end


%calculate power once if we are doing SCV or 1/f method
if sum(strcmp(params.bandMethods,'SCV')+strcmp(params.bandMethods,'PeakPick'))>0
    [~,pow]= multiphasevec2(params.wavefreqs,signal',params.srate,6);
else
    pow = [];
end


%deal with bad data.  instead of removing the time points which makes
%indexing different matrices hard to keep track of, simply replace power values during
%bad times with the NaN and do nanmean when getting bands
%identification
if isfield(params,'bad_data') && ~isempty(pow)
    bad_idx = find(params.bad_data==1);
    pow(:,bad_idx) = NaN;
end


%run the ORCA algorithm using different getband methods and use the best r2
for iL = 1:length(params.bandMethods)
    
    params.getband_method = params.bandMethods{iL};
    
    [PO{iL},fs{iL},bdz{iL},bA{iL},bP{iL},rcn{iL}] = ...
        RunORCA(signal,params,pow);
    
    
    allr2(iL) = rcn{iL}.full_r2;
    %get training and rest test r2
    testr2(iL) = rcn{iL}.crossval_test_r2;
    trainr2(iL) = rcn{iL}.crossval_train_r2;
    
    
    %display r2 for test data
    if isnumeric(params.getband_method) && params.verbose && length(params.getband_method(:))>1 %user defined bands
        disp(['r^2 with user-defined bands:  ',num2str(testr2(iL))]);
    elseif isnumeric(params.getband_method) && params.verbose && length(params.getband_method(:))==1 %random bands
        disp(['r^2 with random bands:  ',num2str(testr2(iL))]);
    
    elseif params.verbose
        disp(['r^2 with ',params.bandMethods{iL},':          ',num2str(testr2(iL))]);
    end
            disp([' '])

end





%determine final outputs after
%we ran everything and we need to extract the best method from r2
bestMethIdx = testr2==max(testr2);

rmfield(params,'getband_method');
NormAmp = PO{bestMethIdx};
bandfreqs = fs{bestMethIdx};
bands = bdz{bestMethIdx};
bandamps = bA{bestMethIdx};
bandphases = bP{bestMethIdx};
recon_stats = rcn{bestMethIdx};

if  isnumeric(params.bandMethods{bestMethIdx}) && length(params.bandMethods{bestMethIdx})==1 %random bands
    recon_stats.best_method = 'Random';

elseif  isnumeric(params.bandMethods{bestMethIdx}) && length(params.bandMethods{bestMethIdx})>1 %user defined bands

    
    recon_stats.best_method = 'User';
else
    recon_stats.best_method = params.bandMethods{bestMethIdx};
end
recon_stats.all_methods =params.bandMethods;
recon_stats.all_r2 = allr2;
recon_stats.testr2=testr2;
recon_stats.trainr2=trainr2;




%extract the r2 for every frequency if we computed
%it using the CoD band id method
CoD_idx = find(strcmp(params.bandMethods,'CoD')==1);
if ~isempty(CoD_idx)
    recon_stats.frequency_r2 = rcn{CoD_idx}.frequency_r2;
    recon_stats.frequency_hz = rcn{CoD_idx}.frequency_hz;
end

%extract the SCV if we computed
%it using the SCV band id method
SCV_idx = find(strcmp(params.bandMethods,'SCV')==1);
if ~isempty(SCV_idx)
    recon_stats.SCV = rcn{SCV_idx}.SCV;
    recon_stats.SCV_freqs = params.wavefreqs;

end


PeakPick_idx = find(strcmp(params.bandMethods,'PeakPick')==1);
if ~isempty(PeakPick_idx)
    recon_stats.PeakPick_pow = rcn{PeakPick_idx}.peakpick_pow;
    recon_stats.PeakPick_freqs = params.wavefreqs;
    recon_stats.PeakPick_fitLine = rcn{PeakPick_idx}.peakpick_fit_line;


end

%determine if we save out all reconstructed traces
% or only the best one
if isfield(params,'saveAllRecons') && params.saveAllRecons
    for iM = 1:length(bA)
        allRecons(iM,:) = rcn{iM}.recon_sig;
        
    end
    recon_stats.recon_sig = allRecons;
end


% output the (perhaps filtered and/or downsampled)
%signal which was analyzed
recon_stats.analyzed_signal = signal;
recon_stats.used_bands = bdz;


%quantify effect size of improvement
if length(allr2)>1
    sortr = sort(fisherstrans(sqrt(allr2)));% sort and convert to z value with fishers transform
    minq = sortr(end)-sortr(end-1);
    maxq = sortr(end)-sortr(1);
    recon_stats.min_CohenQ = minq;
    recon_stats.max_CohenQ = maxq;
else
    recon_stats.min_CohenQ = NaN;
    recon_stats.max_CohenQ = NaN;
end



%display final results and footer regardless of verbose mode
disp(' ')

disp(['Best decomposition with ',recon_stats.best_method,' bands  r^2 = ',num2str(recon_stats.crossval_test_r2)])

if length(allr2)>1
disp(['Effect size (Cohens q): ',num2str(roundn(minq,-2)),' to ',num2str(roundn(maxq,-2))])
end



%footer
disp('************************')
disp(' ')

end

%subfunctions below
%primary ORCA function
function [NormAmp,bandfreqs,bands,bandamps,bandphases,recon_stats] = ...
    RunORCA(signal,params,pow)
wavefreqs = params.wavefreqs;


%step 1, Define or calculate bands. See documentation
%above for all options
if isempty(params.crossval_trainIDX) %empty crossval_training index, so don't do crossval
    if strcmp(params.getband_method,'SCV'); %
        %get SpectCV bands
        [bands,SCV] = Get_SpectCV_Bands(params,pow);
    elseif strcmp(params.getband_method,'CoD');
        [bands,frequency_r2,frequency_hz]  = Get_FreqHoldout_Bands(signal,params);
    elseif isnumeric(params.getband_method) && length(params.getband_method(:))==1 %make random n bands
        n_randbands = (params.getband_method-1);
        bands=[];
        while size(bands,1)<params.getband_method
        bandLims = randi(length(params.wavefreqs),[n_randbands 1]);
        bandLims = unique(sort([min(params.wavefreqs) max(params.wavefreqs) params.wavefreqs(bandLims)]));
        bands(:,1) = bandLims(1:end-1);
        bands(:,2) = bandLims(2:end);
        end
    
    elseif isnumeric(params.getband_method) && length(params.getband_method(:))>1 %user defined bands
        bands = params.getband_method;
        % bands to include full spectrum
        bands(1,1) = min(params.wavefreqs);
        bands(end,2) = max(params.wavefreqs);
    elseif strcmp(params.getband_method,'PeakPick');
       [tmpbands,~,~,~,mean_pow,fit_line] = GetBands(wavefreqs,pow,params.plot_bands);
        bandLims = unique(sort([min(params.wavefreqs) max(params.wavefreqs) tmpbands(:)']));
        bands=[];
        bands(:,1) = bandLims(1:end-1);
        bands(:,2) = bandLims(2:end);     
    end
else %crossvalidate using training indices specified in params
    trainIDX = find(params.crossval_trainIDX==1); 
    testIDX =  find(params.crossval_trainIDX==0);

    
    if strcmp(params.getband_method,'SCV'); %
        %get SpectCV bands
        [bands,SCV] = Get_SpectCV_Bands(params,pow(:,trainIDX));
    elseif strcmp(params.getband_method,'CoD');
        [bands,frequency_r2,frequency_hz]  = Get_FreqHoldout_Bands(signal(trainIDX),params);
     elseif strcmp(params.getband_method,'PeakPick');
        [tmpbands,~,~,~,mean_pow,fit_line] = GetBands(wavefreqs,pow(:,trainIDX),params.plot_bands);

        bandLims = unique(sort([min(params.wavefreqs) max(params.wavefreqs) tmpbands(:)']));
        bands=[];
        bands(:,1) = bandLims(1:end-1);
        bands(:,2) = bandLims(2:end);
        
    elseif isnumeric(params.getband_method) && length(params.getband_method(:))==1 %make random n bands
        n_randbands = (params.getband_method-1);
        bands=[];
        while size(bands,1)<params.getband_method
        
        bandLims = randi(length(params.wavefreqs),[n_randbands 1]);
        bandLims = unique(sort([min(params.wavefreqs) max(params.wavefreqs) params.wavefreqs(bandLims)]));
        
        if (length(bandLims)-1)==params.getband_method
        bands(:,1) = bandLims(1:end-1);
        bands(:,2) = bandLims(2:end);
        
        end
        end
        
    elseif isnumeric(params.getband_method) && length(params.getband_method(:))>1 %user defined bands
        bands = params.getband_method;
        %update bands to include full spectrum
        bands(1,1) = min(params.wavefreqs);
        bands(end,2) = max(params.wavefreqs);
    end
end

%clean up bands and maybe print them
bands = roundn(CheckBands(bands,params.min_bw),-1);
bandStr=num2str(bands);
if params.verbose
    if isnumeric(params.getband_method) && length(params.getband_method(:))>1
       disp(['User bands:'])

        
    elseif isnumeric(params.getband_method) && length(params.getband_method(:))==1 %make random n bands

        disp(['Random bands:'])
    else
        disp([params.getband_method,' bands:'])
    end
    for iB = 1:size(bands,1)
        disp([bandStr(iB,:),' Hz'])
    end
end
  
    %step 2, filter in bands and get amp, phase, and freq (FS)
    bandfreqs = zeros(size(bands,1),length(signal)).*NaN;
    allFilt = zeros(size(bands,1),length(signal)).*NaN;
    bandamps = zeros(size(bands,1),length(signal)).*NaN;
    bandphases = zeros(size(bands,1),length(signal)).*NaN;
    for iB = 1:size(bands,1)
        freq_bands = bands(iB,:);
        [allFilt(iB,:),bandamps(iB,:),bandphases(iB,:),...
            bandfreqs(iB,:)] = FrequencySlide(signal,freq_bands,params);
    end
    
    
    %step 3, Get Normalized Amplitude based on amplitude on a
    %        (half)cycle by cycle analysis
    NormAmp= NaN.*zeros(size(allFilt));
    for iB = 1:size(bands,1)
        [NormAmp(iB,:)] = Get_NormAmp(allFilt(iB,:),bandfreqs(iB,:),params);
    end
    

    %step 4, reconstruct signal to get r2
    [recon_sig,r2] = ReconstructSignal(...
        signal,bands,bandamps,bandphases,bandfreqs);
    
    
    recon_stats.recon_sig = recon_sig;
    recon_stats.full_r2 = r2;
    
    
    %optional step 5, try reconstructing based on shuffled
    %amp phase freq
    if params.numControlPerms
        numShifts = params.numControlPerms;
        for iShift = 1:numShifts
            disp(num2str(iShift))
            shiftAmt = randi([1 size(bandamps,2)]);
            shift_amp = circshift(bandamps,shiftAmt,2);
            shift_phase= circshift(bandphases,shiftAmt,2);
            shift_freq = circshift(bandfreqs,shiftAmt,2);
            
           
            %shift everything an equal amount in time
             [~,shiftAll_r2(iShift)] = ReconstructSignal(...
                 signal,bands,shift_amp,shift_phase,shift_freq);
%         
            %shuffle bands. permute indices of amplitude and phase
            shift_amp = bandamps(randperm(size(bands,1)),:);
            shift_phase= bandphases(randperm(size(bands,1)),:);
            
            %don't shuffle freqs otherwise recon fails==nan
            %shift_freq = bandfreqs(shuffs,:);
            [~,shuffle_r2(iShift)] = ReconstructSignal(...
             signal,bands,shift_amp,shift_phase,bandfreqs);
       
        end
        
        recon_stats.r2_shuffled.time = shiftAll_r2;
        recon_stats.r2_shuffled.bands = shuffle_r2;
  
    else
        recon_stats.r2_shuffled={};
    end
    
        %define r2 on the crossvalidated train and test data
        test_xval_rho = regstats(signal(testIDX),recon_sig(testIDX)','linear','rsquare');
        train_xval_rho = regstats(signal(trainIDX),recon_sig(trainIDX)','linear','rsquare');
        regstats(signal,recon_sig','linear','rsquare');
        recon_stats.crossval_test_r2 = test_xval_rho.rsquare;
        recon_stats.crossval_train_r2 = train_xval_rho.rsquare;

    if strcmp(params.getband_method,'SCV')
        recon_stats.SCV = SCV;
    end   
        
    if strcmp(params.getband_method,'PeakPick')
        
        recon_stats.peakpick_fit_line = fit_line;
        recon_stats.peakpick_pow = mean_pow;
    end   
    
    %get frequency_r2 if we computed it during CoD bandmethod
    if strcmp(params.getband_method,'CoD')
        recon_stats.frequency_r2 = frequency_r2;
        recon_stats.frequency_hz = frequency_hz;
    end   

    %optional save out filtered data
    if params.saveFiltered
        recon_stats.filtered_data = allFilt;
    end

    
end



function [NormAmp] = Get_NormAmp(filtSig,frequencies,params)
%get Normalized amplitude based on percentile of each half cycle
%amplitude

%get distribution of peak and trough amplitudes in filtered signal
[peakampz,peaksamps] = findpeaks(filtSig,1:length(filtSig));
[troughampz,troughsamps] = findpeaks(-filtSig,1:length(filtSig));

%make troughampz back to original sign
troughampz=-troughampz; %need to because the input was flipped sign above


%declare first datapoint as either peak or trough
if peaksamps(1)<troughsamps(1)
    troughsamps = [1 troughsamps];
    troughampz = [filtSig(1) troughampz];
else
    peaksamps = [1 peaksamps];
    peakampz = [filtSig(1) peakampz];
end

%     %declare last datpoint as either peak or trough
if peaksamps(end)>troughsamps(end)
    troughsamps = [troughsamps length(filtSig)];
    troughampz = [troughampz filtSig(end) ];
else
    peaksamps = [peaksamps length(filtSig)];
    peakampz = [peakampz filtSig(end)];
end

allsamps = sort([peaksamps troughsamps]);
allVoltages = filtSig(allsamps);
for iS = 1:length(allsamps)-1
    waveAmp(allsamps(iS):allsamps(iS+1)) = ...
        abs(allVoltages(iS)-allVoltages(iS+1));
end

% default:normalize amplitude by instantaneous frequency 
% This will ensure that activity at slower frequencies
% in the band are not always regarded as larger due to
% the general 1/f nature of amplitude.
if params.normalize_amp
tmpNA = waveAmp.*frequencies;
else
    tmpNA = waveAmp; %don't normalize.
end

NormAmp = GetPRCT(tmpNA);
NormAmp(isnan(NormAmp)) = 0; %replace NaNs with zero
end

function [filtered_signal,bandamps,bandPhase,bandFreq] = FrequencySlide(signal,freq_bands,params)


%filter using variable transition width depending on
%frequency band. Try several and keep final output as that which
%is most similar to raw signal
transWidths = .01:.03:.15;
for iI = 1:length(transWidths)
    trans_width = transWidths(iI);
    idealresponse  = [ 0 0 1 1 0 0 ];
    filtfreqbounds = [ 0 (1-trans_width)*freq_bands(1) freq_bands(1) freq_bands(2) freq_bands(2)*(1+trans_width) params.srate/2 ]/(params.srate/2);
    filt_order     = round(2*(params.srate/freq_bands(1))); %modified this to deal with weirdness when this number got above ~ 6000samples
    filterweights  = firls(filt_order,filtfreqbounds,idealresponse);
    tmp_filtered_signal(iI,:) = filtfilt(filterweights,1,signal);
    thisCorr(iI) = corr(signal,tmp_filtered_signal(iI,:)');
end

%keep the filtered signal with the highest correlation
%to the raw signal. tie goes to largest trans width
best_transwidth = max(find(thisCorr==max(thisCorr)));
filtered_signal = tmp_filtered_signal(best_transwidth,:);


%hilbert the filtered signal
temphilbert = hilbert(filtered_signal);
anglehilbert = angle(temphilbert);

bandPhase = anglehilbert;
bandamps = abs(temphilbert);

%code from frequency sliding algorithm (Cohen, 2014 J Neuroscience)
frompaper = params.srate*diff(unwrap(anglehilbert))/(2*pi); 
frompaper(end+1) = NaN; %deal with fact that diff loses a sample
time_wins = [.05 .2 .4]; %time windows in fractions of a second
orders = [time_wins*params.srate];
% %median filter using different window sizes. window signal to make it more
% %tractable
numchunks = 10;
chunks = floor(linspace(1,length(frompaper),numchunks));
meds = zeros(length(orders),length(frompaper));
for iWin = 1:length(orders)
    for iChunk = 2:numchunks
        chunkidx = chunks(iChunk-1):chunks(iChunk)-1;
        foo = medfilt1(frompaper(chunkidx),round(orders(iWin)));      
        meds(iWin,chunkidx) = foo;
    end
end

%take the median value, cohen method
median_of_meds = median(meds);

%NaN out values outside of the filter band

clear below* above* outside*
below_idx = (median_of_meds<freq_bands(1));
above_idx = (median_of_meds>freq_bands(2));
outside_idx = find([below_idx+above_idx]==1);
median_of_meds(outside_idx)=NaN;
bandFreq = median_of_meds; %final frequency sliding estimate

end

function [percentile] = GetPRCT(observations)
%gets percentile of all observations in a distribution.
%input: observations is a 1Xn vector
%output: percentile is the percentile for each value in observations

%get percentile of all non nan values
non_nan = sum(~isnan(observations));
percentile = tiedrank(observations)/non_nan;
end

function  [freq_bands,bandidx,bandpow,fs_nantimes,mean_pow,fit_line] = GetBands(wavefreqs,pow,plotfit)
%This is the original formulation of "peak picking" or 1/f slop fitting
%described in Watrous et al., 2018 eLife

%fs nantimes are the times where average band power is below 1/f fit line
%in each band
fz = log(wavefreqs);
mean_pow = log(nanmean(pow,2)); %nanmean nov 1 to deal with pow values that are nan based on bad data
[b,~] = robustfit(fz,mean_pow);
fit_line = b(1)+b(2).*fz;

above1f = (mean_pow-fit_line')>0;
bw = bwlabel(above1f);
ctr=1;

for iB = 1:max(unique(bw))
    %make sure its actually a band and not a point frequency
    idx = find(bw==iB);
    if length(idx)>1
        freq_bands(ctr,1) = wavefreqs(min(idx));
        freq_bands(ctr,2) = wavefreqs(max(idx));
        bandidx{ctr} = idx;
        bandpow(ctr,:) = log(mean(pow(idx,:)));
        crit_pow = mean(fit_line(idx));
        fs_nantimes{ctr} =find(bandpow<crit_pow);
        ctr=ctr+1;
    end
end


if plotfit
    
    plotx = [1 10:20:max(wavefreqs)];
    figure('Position',[100 100 1200 600]);
    %subplot(8,1,2:5)
    plot(fz,mean_pow)
    hold on;
    plot(fz,fit_line,'--r')
    fooo = dsearchn(wavefreqs',plotx');
    set(gca,'YTick',[],'XTick',fz(fooo),'XTickLabel',plotx);
    xlabel('Frequency (Hz)')
    ylabel('Mean Log Power')
    % title('Global Fit')
    axis tight
    xlim(log([min(wavefreqs) max(wavefreqs)]))
    yl =ylim;
    colr = linspecer(size(bandidx,2));
    for iB = 1:size(bandidx,2)
        xx = fz([min(bandidx{iB}) max(bandidx{iB}) max(bandidx{iB}) min(bandidx{iB})]);;
        yy = [yl(1) yl(1) yl(2) yl(2)];
        patch(xx,yy,colr(iB,:),'FaceAlpha',.3)
    end
    legend({'Signal','Fit Line'})
    title('Bands defined by PeakPick method')
    pubgraph_AW({'fs','lw'},{14,2})
    
end

end

function [bands,spectCV] = Get_SpectCV_Bands(params,pow);
%calculate spectral coefficient of variation and use local minima in this 
%function (across frequencies) to define bands.

%calculate spectral coefficient of variation
spectCV = std(pow')./mean(pow,2)';
iSpectCV = spectCV.*-1; %invert it so find peaks actually finds troughs

%find minima of spectral COV function
[~,minlocs] = findpeaks(iSpectCV,params.wavefreqs,'MinPeakDistance',params.min_bw);
minlocs = [min(params.wavefreqs) minlocs]; %set lowest frequency as a boundary

ctr=1;
for iMin =1:length(minlocs)
    if iMin<length(minlocs)
        f1 = minlocs(iMin);
        f2 = minlocs(iMin+1);
        
    else
        f1 = minlocs(iMin);
        f2 = max(params.wavefreqs);
    end
    bands(ctr,1) = f1;
    bands(ctr,2) = f2;
    ctr=ctr+1;
end


%ensure bands are good
bands = CheckBands(bands,params.min_bw);

if isfield(params,'plot_bands') && params.plot_bands %plot bands?
    
    figure('Position',[100 100 1200 600]);
    [yy,h1,h2]=plotyy(log(params.wavefreqs),log(mean(pow,2)),log(params.wavefreqs),spectCV);
    ylabel(yy(1),'Log Mean Power');ylabel(yy(2),'Coefficient of Variation');
    yl =ylim;
    xlim(yy,log([min(params.wavefreqs) max(params.wavefreqs)]))
    
    colr = linspecer(size(bands,1));
    for iB = 1:size(bands,1)
        xx = log([bands(iB,1) bands(iB,2) bands(iB,2) bands(iB,1)]);;
        yy = [yl(1) yl(1) yl(2) yl(2)];
        patch(xx,yy,colr(iB,:),'FaceAlpha',.3)
    end
    xt=get(gca,'XTick');
    set(gca,'XTickLabel',roundn(exp(xt),-1))
    title('Bands defined by SCV method')  
    pubgraph_AW({'fs','lw'},{14,2})
end
end

function [bands,frequency_r2,frequency_hz] = Get_FreqHoldout_Bands(signal,params);
% get band edges by looking at % explained variance in signal
% reconstruction when including vs excluding individual frequencies
% bands are defined based on local minima in "frequency_r2" along with
% placing a band edge where frequency_r2<(1/length(frequency_hz)
% use params.plot_bands=1 to see how this works in more detail

%outputs bands and the explained variance of activity
%at each frequency (frequency_r2) as well as the tested frequencies
%(frequency_hz)

%do cwt
[wt,frequency_hz] = cwt(signal,params.srate);

%crop to frequencies in range
f_idx = intersect(find(frequency_hz>=min(params.wavefreqs)),...
    find(frequency_hz<=max(params.wavefreqs))); %define indices of frequencies to keep
wt = wt(f_idx,:);
frequency_hz = frequency_hz(f_idx);

%flip it ( even though this doesn't impact results based on how icwt.m works)
frequency_hz = flipud(frequency_hz);
wt=flipud(wt);

%compute r2 for each frequency by zeroing activity at other freqs,
%reconstructing signal, and computing fit to original signal
frequency_r2 = ones(1,length(frequency_hz)).*NaN;
for iFr = 1:length(frequency_hz)
    omit_wt = wt;
    omit_idx = setdiff(1:length(frequency_hz),iFr);
    omit_wt(omit_idx,:) = 0;
    recon_sig = icwt(omit_wt);
    frequency_r2(iFr) = corr(signal,recon_sig')^2;
end


%define bands based on local minima of explained variance
[~,idx2] = findpeaks(-frequency_r2,frequency_hz,'MinPeakDistance',params.min_bw);

%define band boundary as place where explained variance drops below what is
%expected if each individual contributes equally (i.e. 1/length(frequency_hz)
critical_r2 = 1/length(frequency_hz);
bwww = bwlabel(frequency_r2<critical_r2);
idx3 = [];
for iClus = 1:max(bwww)
    idx3 = [idx3 frequency_hz(min(find(bwww==iClus)))];
end

minlocs = sort([min(params.wavefreqs); idx2; idx3']); %include lowest freq as a boundary
ctr=1;
for iMin =1:length(minlocs)
    if iMin<length(minlocs)
        f1 = minlocs(iMin);
        f2 = minlocs(iMin+1);
        
    else
        f1 = minlocs(iMin);
        f2 = max(params.wavefreqs);
    end
    bands(ctr,1) = f1;
    bands(ctr,2) = f2;
    ctr=ctr+1;
end

%ensure bands are good
bands = CheckBands(bands,params.min_bw);

if isfield(params,'plot_bands') && params.plot_bands %maybe plot
    
    figure('Position',[100 100 1200 600]);
    plot(log(frequency_hz),frequency_r2)
    ylabel('Explained Variance (r^2)');
    yl =ylim;
    
    colr = linspecer(size(bands,1));
    for iB = 1:size(bands,1)
        xx = log([bands(iB,1) bands(iB,2) bands(iB,2) bands(iB,1)]);;
        yy = [yl(1) yl(1) yl(2) yl(2)];
        patch(xx,yy,colr(iB,:),'FaceAlpha',.3)
    end
    xlim(log([min(params.wavefreqs) max(params.wavefreqs)]))
    xlabel('Frequency (Hz)')
    xt=get(gca,'XTick');
    set(gca,'XTickLabel',roundn(exp(xt),-1));
    title('Bands defined by CoD method')
    axis tight
    pubgraph_AW({'fs','lw'},{14,2})
end

end

function bands = CheckBands(inbands,min_bw)
%only keep bands that are wider than min_bw
%this fixes issue with too narrow bands leading to
%failure to filter
bands = inbands(diff(inbands')>min_bw,:);
end


%reconstruct signal using band amp phase freq
function [recon_sig,r2] = ReconstructSignal(sig,bands,bandamps,bandphases,bandfreqs)
%outputs the reconstructed signal,
%the r2 fit to the real signal,

%round frequency sliding to make it tractable &
%filter out frequencies outside of our bands. Rounding to more
%precise frequencies e.g. roundn(freqs,-1) increases
%r2 very marginally (test case went from .9271 to .9228)
%but increases computation time drastically (376s to 11s)

id1 = bandfreqs<min(bands(:));
id2 = bandfreqs>max(bands(:));
thresh_IF = roundn(bandfreqs,-1); %round it to nearest .1 Hz to speed computation
%thresh_IF = round(bandfreqs); %round it to nearest 1 Hz
thresh_IF(id1) = NaN;thresh_IF(id2)= NaN;

unFz = unique(thresh_IF);
unFz = unFz(~isnan(unFz)); %crop NaNs

%generate logspaced frequency indices.  generate many so we
%have a high chance of capturing all of the unique instances
%of values in bandfreqs
fz = logspace(log10(min(bands(:))),log10(max(bands(:))),length(unFz).*1000)'; %transpose so column dim ==1 for dsearchn

%ensure its N X 1
if size(unFz,2)>1
    unFz=unFz';
end
try
    check_fz_idx = dsearchn(fz,unFz); %both should be N X 1
catch
end

unique_found_freqs = unique(check_fz_idx);
%check to make sure there are equal matches between fz and unFz before
%continuing
if length(unique_found_freqs)<length(unFz);
    warning('Failed to match frequencies')
end

%generate final set of frequencies to be used
fz = fz(unique_found_freqs);

%make pseudo cwt.m output matrix using amplitude and phase
%of detected activity
complexMat = zeros(length(fz),length(sig));

%build complexMat
for iFr = 1:length(fz)
    
   %which band is this freq in?
    band_idx = find((fz(iFr)>=bands(:,1))+(fz(iFr)<=bands(:,2))==2);
   
    %which samples have activity at this frequency?
    foundFreqIDX =   find(thresh_IF(band_idx,:)==roundn(fz(iFr),-1));

    %make complex values from amplitude and phase
    c =bandamps(band_idx,foundFreqIDX).*exp(bandphases(band_idx,foundFreqIDX)*1i);
    if all(~isinf(c)) %if all values aren't inf
        complexMat(iFr,foundFreqIDX) = c; %put the values in the right spot
    else %replace inf values with 0, in very rare cases
        nInf = sum(isinf(c));
        c(isinf(c))=0;
        disp(['Found and replaced ',num2str(nInf),' Inf with 0 during reconstruction'])       
    end
end

%reconstruct signal and calculate R2
recon_sig = icwt(complexMat);
r2=regstats(sig,recon_sig','linear','rsquare') ;
r2 = r2.rsquare;

end

