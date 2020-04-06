---
title: "Frequency analysis"
author: "Lau Møller Andersen & Mikkel C. Vinding"
date: "Feb. 2018; NatMEG, Stockholm, Sweden"
output: html_document
---

### Set up general paths
Change these to the appropriate paths for your operating system and setup.

```{r, engine='octave', eval=FALSE}
addpath /home/lau/matlab/fieldtrip-20170613_workshop_natmeg/
ft_defaults

raw_meg_path = '/archive/20067_workshop_source_reconstruction/MEG/';
meg_path = '/home/share/workshop_source_reconstruction/data/MEG/';
mri_path = '/home/share/workshop_source_reconstruction/data/MRI/';

%% subjects and dates
subjects_and_dates = ...
                    {
                        'NatMEG_0177/170424/'
                    };
        
data_path = fullfile(meg_path, subjects_and_dates{1});
cd(data_path);
```

### Load data
Read in cleaned data from yesterday (remember where you put the data and what you named it).

If you did not complete the data preperation tutorial, you can load the data file _cleaned_downsampled_data.mat_ from the tutorial material:

```{r, engine='octave', eval=FALSE}
load('cleaned_downsampled_data.mat');
disp('Done');

```

In the following tutorial, we will only analyse the conditions where the index finger was stimulated, i.e. the condition corresponding to trigger value 8. Use _ft_selectdata_ to select the conditions with trigger 8:

```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.trials = cleaned_downsampled_data.trialinfo==8;

epochs = ft_selectdata(cfg, cleaned_downsampled_data)
```


First, we baseline correct the data that we are going to use for the TFRs - i.e. demean the epochs.

```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.demean          = 'yes';
cfg.baselinewindow  = 'all';

epochs_bs = ft_preprocessing(cfg, epochs);

```

## Estimate Power Spectral Density
We are going to compute the power spectral density (PSD) across the epochs. This will give us the average power across all trials.
We will use different methods to estimate PDS: 1) a single taper window, and 2) multi-tapered windows. All calculations will be done with the FieldTrip function _ft_freqanalysis_. ft_freqanalysis is the default function for running various types of frequency analysis such as calculating PSD and doing the time-frequency analysis. To get more info type:


```{r, engine='octave', eval=FALSE}
help ft_freqanalysis
```


###  Get PSD with single taper
In the first calculation, we will use a Hann window to taper the epochs. This is applied to each epoch in the data structure before doing Fourier decomposition and then calculating the power. The Hann window is selected by passing the configuration _cfg.taper = 'hanning'_ to ft_freqanalysis.

NB. If the process is slow, change _cfg.channel_ to 'MEG*1' to select only MEG magnetometers to run. You can always come back and re-run the analysis.


```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.output      = 'pow';          % Return PSD
cfg.channel     = {'MEG','EEG'};  % Calculate for MEG and EEG
cfg.method      = 'mtmfft';
cfg.taper       = 'hanning';      % Hann window as taper
cfg.foilim      = [1 95];         % Frequency range

psd_hann = ft_freqanalysis(cfg, epochs_bs);

```

Once finished, look what is in the structure _psd_hann_.

* Q: What is the dimension of the data, and what do the dimensions represent?

Plot the PSD using _ft_multiplotER_. the configuration _cfg.xlim_ specifies the limits of the frequency axis.


```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = 'neuromag306mag'; % Layout for MEG magnetometers
cfg.showlabels      = 'yes';
cfg.xlim            = [3 35];           % Frequencies to plot

figure;
ft_multiplotER(cfg, psd_hann);

```

![](./images/PSD_hann.jpg)

### Get PSD with multitapers
Now we will do the same analysis, but use tapers based on Slepian sequences of tapers. The spectral smoothing is specified by passing the frequency range in the configuration _cfg.tapsmofrq_ when calling _ft_freqanalysis_. From this _ft_freqanalysis_ will calculate the amount of tapers needed to achieve the desired frequency smoothing. Note that the frequency smoothing is the range in either direction, so the frequency smoothing is the double of the numerical value given in _cfg.tapsmofrq.

* Q: How many tapers is used when you run the command below? Hint: it tells you in the terminal.


```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.output      = 'pow';          % Return PSD
cfg.channel     = {'MEG','EEG'};  % Calculate for MEG and EEG
cfg.method      = 'mtmfft';
cfg.taper       = 'dpss';         % Multitapers based on Slepian sequences
cfg.tapsmofrq   = 2;              % Smoothing +/- 2 Hz
cfg.foilim      = [1 95];

psd_dpss = ft_freqanalysis(cfg, epochs_bs);

```

Plot the PSD using _ft_multiplotER_ as above.

Try to change the frequency smoothing range to e.g 10 Hz.

```{r, engine='octave', eval=FALSE}
cfg.tapsmofrq   = 10;            % Q: How many tapers does this use?
psd_dpss10 = ft_freqanalysis(cfg, epochs_bs);

```

Compare the single taper PSD, and the two multitaper PSD you have calculated. Plot them side-by-side using _ft_multiplotER_

```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = 'neuromag306mag'; % Layout for MEG magnetometers
cfg.showlabels      = 'yes';
cfg.xlim            = [3 35];           % Frequencies to plot

figure;
ft_multiplotER(cfg, psd_hann, psd_dpss, psd_dpss10);

```

![](./images/PSD_all.jpg)

Compare the results from the different methods to calculate PSD: Select alpha range (~8-12 Hz) in the multiplot to plot as topoplots. 

* Q: How different/alike are they? Explain in your own words why?

Try to select the beta range (~14-30 Hz) and compare topoplots.

* Q: How different/alike are they? Explain in your own words why?

Try to plot the "high-gamma" range (~55-95 Hz) by changing _cfg.xlim_ to [55 95] and call ft_multiplotER_ as before.

* Q: How do the high-gamma spectra compare between methods? Explain in your own words why?

Bonus: Change the layout to view the PSD of EEG

```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.parameter   = 'powspctrm';
cfg.layout      = 'neuromag306eeg1005_natmeg.lay'; # EEG layout
cfg.showlabels  = 'yes';
cfg.xlim        = [3 35];
figure;
ft_multiplotER(cfg, psd_hann);
```



## Time-frequency analysis
The PSD analysis above (sort of) assumes that the spectral power is the same across the entire epoch. This is the way we often analyze resting state data or similar paradigms. But we can also analyze how the spectral signals evolve over time, e.g. how oscillatory activity changes as a response to stimuli.

Instead of calculating the power across the entire epoch and then average all epochs, we now will calculate the power for each time sample: We center the window on a given time sample and estimate the power of that particular window. As above, we can choose to taper the window with a single taper (e.g., the Hann window) or by using multitapers. Finally, we will also use Wavelet analysis, which are sequences of tapered sine-waves fitted to the data for the given time point.

### Get TFR with single taper MEG
In FieldTrip we calculate the time-frequency response (TFR) with the function _ft_freqanalysis_, as we did for the PSD. This time we change the method to 'mtmconvol'. We also have to specify the time sample resolution (_toi_ -- times of interest) and the length of the window centered on the time of interest (_t_ftimwin_ -- window length).

NB. TFR calculations take significantly longer time to calculate than PSD. You can always change _cfg.channel_ to  'MEG*1' to save time. You can always go back and redo the other channel types (but do make sure you finish the exercise).


```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.output      = 'pow';      
cfg.channel     = {'MEG','EEG'};      % Change to 'MEG*1' to analyse all channels
cfg.method      = 'mtmconvol';
cfg.taper       = 'hanning';    % Hann window as taper
cfg.foi         = 1:1:45;       % Frequencies we want to estimate from 1 Hz to 45 Hz in steps of 1HZ
cfg.toi         = -1.650:0.01:1.650;            % Timepoints to center on
cfg.t_ftimwin   = 0.5*ones(length(cfg.foi),1);  % length of time window

tfr_hann = ft_freqanalysis(cfg, epochs_bs);
```

Plot the result with _ft_multiplotTFR_:

```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = 'neuromag306mag';
cfg.showlabels      = 'yes';

figure;
ft_multiplotTFR(cfg, tfr_hann);
```

* Q: This plot look weird! How come? What can we see from the plot (Hint: Remember the PSD plots from before)?

![](./images/TFR_noBaseline.jpg)

Do the same plot, but this time define a baseline to use as a reference for plotting. This can be from the start of the trial to stimuli at time 0 [-inf 0], or the entire epoch [-inf inf]. You can see the different options by calling _help_ for _ft_multiplotTFR_. Here we plot the relative change (_cfg.baselinetype = 'relative'_) compared to the baseline from the start to time 0.


```{r, engine='octave', eval=FALSE}
% Plot
cfg = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = 'neuromag306mag';
cfg.showlabels      = 'yes';
cfg.baselinetype    = 'relative';  % Type of baseline, see help ft_multiplotTFR
cfg.baseline        = [-inf 0];    % Time of baseline

figure; ft_multiplotTFR(cfg, tfr_hann);
````

![](./images/TFR_hannBaseline.jpg)

Zoom in on a single channel:

![](./images/TFR_hannSingleChan.jpg)

* Q: How does the baseline change what we can infer from the TFR (Hint: Toggle the colorbar option to show what the colors represent in either plot)?

### Get TFR with Hann taper with varying length
The previous TFR analysis used a fixed window for all frequencies. This is not optimal in this case. Now we will run the same analysis, but with time windows that vary with the individual frequency. Here we use a time window that corresponds to five times the wavelength of frequency of interest (_cfg.t_ftimwin = 5./cfg.foi_)

Before proceeding, compare the numeric values of _cfg.t_ftimwin_ from the previous to the current.

```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.output      = 'pow';        % Return PSD
cfg.channel     = {'MEG','EEG'};      % Change to 'MEG*1' to analyse all channels
cfg.method      = 'mtmconvol';
cfg.taper       = 'hanning';    % Hann window as taper
cfg.foi         = 1:1:45;       % Frequencies we want to estimate from 1 Hz to 45 Hz in steps of 1HZ
cfg.toi         = -1.650:0.01:1.650;             % Times to center on
cfg.t_ftimwin    = 5./cfg.foi; 

tfr_hann5 = ft_freqanalysis(cfg, epochs_bs);
````

Plot the results as before (remember the baseline):

```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = 'neuromag306mag';
cfg.showlabels      = 'yes';
cfg.baselinetype    = 'relative';
cfg.baseline        = [-inf 0];

figure; ft_multiplotTFR(cfg, tfr_hann5);
````

![](./images/TFR_hann5.jpg)

Zoom in on a single channel:

![](./images/TFR_hann5singChan.jpg)


* Q: Is something missing---why the round edges of the plot?


### Get TFR with multitaper
Now we will do the same with multitapers. The smoothing range is given in _cfg.tapsmofrq_. Take a look at the values _cfg.tapsmofrq_ and _cfg.t_ftimwin_ before running _ft_freqanalysis_. 

Based on the values, how do you expect the plot of the results to look?

```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.output      = 'pow';        
cfg.channel     = {'MEG','EEG'};      % Change to 'MEG*1' to analyse all channels
cfg.method      = 'mtmconvol';
cfg.taper       = 'dpss';       % Slepian sequence as tapers
cfg.foi         = 1:1:45;       % Frequencies we want to estimate from 1 Hz to 45 Hz in steps of 1HZ
cfg.toi         = -1.650:0.01:1.650;             % Times to center on
cfg.t_ftimwin   = 5./cfg.foi;   % length of time window
cfg.tapsmofrq   = 0.5 *cfg.foi; % Smoothing

tfr_dpss = ft_freqanalysis(cfg, epochs_bs);
````

Plot the results as before with _ft_multiplotTFR_. Remember the baseline.

Zoom in on single channels. Do the plot compare to your expectation?

![](./images/TFR_dpss.jpg)


### TFR with Morlet wavelets
Wavelets are different from the previous, in that the "tapering" is done on a wave-function that is fitted to the signal centered on the time-point of interest---in contrast to applying the taper on the data itself as you did before.

To do wavelet analysis in FieldTrip, change the _cfg.method_ configuration passed to _fr_freqanalysis_ to 'wavelet'. Again, we specify the timepoint of interest (_cfg.toi_), where we center the wavelets and the frequency of interest (_cfg.foi_) as before. For wavelets, we also need to specify the "width" of the wavelet functions, i.e., the number of cycles in each function.


```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.channel     = 'MEG*1';
cfg.method      = 'wavelet';
cfg.foi         = 1:1:45;       % Frequencies we want to estimate from 1 Hz to 45 Hz in steps of 1HZ
cfg.toi         = -1.650:0.01:1.650;
cfg.width       = 5;
cfg.pad         = 'nextpow2';

tfr_wavelet = ft_freqanalysis(cfg, epochs_bs);
````

Plot the results with _ft_multiplotTFR_ as before.

![](./images/TFR_wav.jpg)


Plot the EEG TFRs for comparison:
```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = 'neuromag306eeg1005_natmeg.lay';
cfg.showlabels      = 'yes';
cfg.baselinetype    = 'relative';
cfg.baseline        = [-inf inf];

figure; ft_multiplotTFR(cfg, tfr_wavelet);
````

### Increase sensitivity to non-phase locked signals
One of the purposes of TFR analysis is to reveal signals that are non-phase locked to stimuli, in contrast to the analysis of evoked fields/potentials. But the phase-locked signals also show up in the TFRs, as they are present in the signal and therefore decomposed by the Fourier transform. 
To enhance the sensitivity to the non-phase locked activity, we can subtract the phase locked activity and redo the analysis.

Before you proceed, guess how we can find the phase locked part of the signal (Hint: you might already have done that)

The phase-locked activity is the evoked field/potential. Assuming the average evoked signal are a good representation of the phase-locked signal (i.e., sufficient signal-to-noise ratio), we can subtract it from each trial, and redo the TFR analysis.


```{r, engine='octave', eval=FALSE}
% Get timelocked again
timelock = ft_timelockanalysis([],epochs_bs);

epochs2 = epochs_bs;
for i = 1:length(epochs_bs)
    epochs2.trial{i} = epochs_bs.trial{i} - timelock.avg;
end

% TFR with Morlet wavelets
cfg = [];
cfg.method      = 'wavelet';
cfg.foi         = 1:1:45;       % Frequencies we want to estimate from 1 Hz to 45 Hz in steps of 1HZ
cfg.toi         = -1.650:0.01:1.650;
cfg.width       = 7;
cfg.pad         = 'nextpow2';

tfr_wavelet2 = ft_freqanalysis(cfg, epochs2);
```

Plot the new TFR. Can you spot the difference from the previous wavelet analysis?

```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = 'neuromag306mag';
cfg.baselinetype    = 'relative';
cfg.baseline        = [-inf 0];
cfg.zlim            = [.5 1.8]        % Make sure the plots have same color scaling

figure; ft_multiplotTFR(cfg, tfr_wavelet);   % First wavelet analysis
figure; ft_multiplotTFR(cfg, tfr_wavelet2);  % Second wavelet analysis
```

![](./images/TFR_wav2.jpg)



### Final question
Now you have compared different methods to calculate TFRs.

* Q: What do the results show? Pick the method you prefer and explain what type of induced responses we see (Hint: Remember what ERS and ERD stood for). It is a good idea to use representative plots to illustrate the results .


