

In this tutorial, we will do source reconstruction of time-frequency data with a beamformer method known as Dynamic Imaging of Coherent Sources (DICS; Gross et al. 2001). DICS is based on reconstructing sources that show strong dependency (coherence) in the frequency domain.

We will start by finding a time-frequency area of interest that want to know the underlying sources. From time-domain data, we will calculate the cross-spectral density of the frequency (or frequencies) of interest and then use the cross-spectral density to find underlying the sources. Most of this is done "under the hood" by FieldTrip with the function _ft_sourceanalysis_.

First, we need to prepare the raw data, define a source model, and calculate the lead field.

### Set up general paths
Note that if work from the path where all the downloadable parts are downloaded to, you don't need to change the paths(leave them as is, but do evaluate the sections). The paths serve as an example of how you can set up your analysis structure, which is especially useful if you have more than one subject 

Change these to appropriate paths for your operating system and setup. Note that if work from the path where all the downloadable parts are downloaded to, you don't need to change the paths (leave them as is, but do evaluate the section).

```octave
addpath /home/lau/matlab/fieldtrip-20170613_workshop_natmeg/
ft_defaults

raw_meg_path = '/archive/20067_workshop_source_reconstruction/MEG/';
meg_path = '/home/share/workshop_source_reconstruction/data/MEG/';
mri_path = '/home/share/workshop_source_reconstruction/data/MRI/';

```

Make subject and recording specific paths (the cell array "subjects_and_dates" can be expanded)

```{r, engine='octave', eval=FALSE}
%% subjects and dates

subjects_and_dates = ...
                    {
                        'NatMEG_0177/170424/'
                    };

output_path = fullfile(meg_path, subjects_and_dates{1});

events = [1 2 4 8 16]; %% right little, right ring, right middle, right index, right thumb

```

### Load files necessary for the beamforming

We are loading four different files here:  

1. We are loading the TFRs to locate interesting targets for beamformer source reconstruction  
2. We are loading the head model for the MEG data  
3. We are loading the head model for the EEG data  
4. We are loading the preprocessed data (baseline_data) on which we will do the actual source reconstruction  

If you saved the files from the preprocessing tutorial, head model tutorial, and TFR analysis tutorial you should load those files. If not, you can find the following files in the tutorial material folder:

```{r, engine='octave', eval=FALSE}
%% go to relevant path and load data

cd(output_path)
disp 'Loading input data'
load combined_tfrs.mat
load headmodel_meg.mat
load headmodel_eeg.mat
load baseline_data.mat
load cleaned_downsampled_data.mat
disp Done
```

### Set channels
Here you can set which channels you want to base the source reconstruction on. For now, we will use the gradiometers, but you can always go back and try to re-run the tutoral with the other sensor types. You only need to change <channels> to one of the three options specified  

```{r, engine='octave', eval=FALSE}
%% set channels

channels = 'meggrad';  % either 'megmag', 'meggrad' or 'eeg'

if strcmp(channels, 'eeg')
    headmodel = headmodel_eeg;
    sensor_type = 'eeg';
else
    headmodel = headmodel_meg;
    sensor_type = 'meg';
end
```

### Identify interesting features

Plot all the channels in the data. Notice (at least) three features of interest. Please explore the plots!  

1. The so-called beta rebound from about 700 to 1000 msec at around 16 Hz  
2. The so-called mu desynchronization from about 300 to 800 msec at around 10 Hz  
3. The so-called high-beta desynchronization from about  200 to 500 msec at around 22 Hz

```{r, engine='octave', eval=FALSE}
%% identify interesting features in the data

% pick an event
event_index = 1;

cfg = [];
cfg.layout = 'neuromag306cmb.lay'; %% this is combined gradiometers, you can also choose <'neuromag306mag.lay'> or <'neuromag306eeg1005_natmeg.lay'
cfg.baselinetype = 'relative'; %% absolute power is not terribly meaningful; therefore we use 'relative' to look at increases and decreases in power relative to the overall power level
cfg.baseline = [-Inf Inf]; %% from min to max
cfg.colorbar = 'yes'; %% show the interpretation of the colours
% cfg.zlim = [0.6 1.6]; %% play around with this parameter to familiarize yourself with these plots

ft_multiplotTFR(cfg, combined_tfrs{event_index});
```


![](./images/tfr_example_tactile.png)

![](./images/tfr_example_occipital.png)

### Visualizing channels
```{r, engine='octave', eval=FALSE}
%% channels that we'll plot throughout

colours = ones(306, 3);
tactile_channel = 'MEG0432';
occipital_channel = 'MEG2142';
tactile_channel_index = find(strcmp(baseline_data.label, tactile_channel));
occipital_channel_index = find(strcmp(baseline_data.label, occipital_channel));

colours(tactile_channel_index, :)   = [1 0 0]; %% make red
colours(tactile_channel_index + 1, :)   = [1 0 0]; %% make red
colours(occipital_channel_index, :) = [0 0 1]; %% make blue
colours(occipital_channel_index + 1, :) = [0 0 1]; %% make blue

figure('units', 'normalized', 'outerposition', [0 0 1 1]);
hold on
sensors = cleaned_downsampled_data.grad;
ft_plot_sens(sensors, 'facecolor', colours, 'facealpha', 0.8);

ft_plot_vol(ft_convert_units(headmodel_eeg, 'cm'));

view([-45 25])
```

### Showing where the channels are (we will use these throughout)
![](./images/sensors_and_brain.png)

## Finding the source of the "beta rebound"
We will here focus on reconstructing the activity underlying the beta rebound.  

1. We crop the data to find the period of interest (690 to 970 msec).
2. We define a baseline period of similar duration. We will compare the period of interest against this since the power estimates that the beamformer results in are not informative in themselves.

In this step, we will also make a combination of the time-window of interest and the baseline. Its use will become apparent later.

```{r, engine='octave', eval=FALSE}
%% time window of interest Beta Rebound

beta_toi     = [0.690 0.970];
baseline_toi = [-0.500 -0.220];
n_events = length(events);

tois_rebound = cell(1, n_events);
tois_baseline  = cell(1, n_events);

% Run for all events
for event_index = 1:n_events
    
    event = events(event_index);
    
    cfg = [];
    cfg.toilim = beta_toi;
    cfg.trials = baseline_data.trialinfo == event;
    
    tois_rebound{event_index} = ft_redefinetrial(cfg, baseline_data);
    
    cfg.toilim = baseline_toi;
    
    tois_baseline{event_index} = ft_redefinetrial(cfg, baseline_data);
    
end

% combined data
tois_combined = cell(1, n_events);

for event_index = 1:n_events

    cfg = [];

    tois_combined{event_index} = ft_appenddata(cfg, tois_rebound{event_index}, tois_baseline{event_index});
    
end
```

### Example plot time-locked

To get an idea about the new data snippets, let us take a look at the data. The main lesson here is that there is no time-locked activity.

```{r, engine='octave', eval=FALSE}
%% example plot tois

event_index = 1;
channel_of_interest = 'MEG0431';
toi_rebound  = tois_rebound{event_index};
toi_baseline = tois_baseline{event_index};
channel_index = strcmp(channel_of_interest, toi_rebound.label);
n_samples = length(toi_rebound.time{1});
n_trials = length(toi_rebound.trial);

figure('units', 'normalized', 'outerposition', [0 0 1 1]);
hold on
for n_trial = 1:n_trials
    plot(1:n_samples, toi_rebound.trial{n_trial}(channel_index, :), 'r')    
end
xlabel('Sample no');
ylabel('Magnetic Field Strength (T)');
title(channel_of_interest)
for n_trial = 1:n_trials
    plot(1:n_samples, toi_baseline.trial{n_trial}(channel_index, :), 'b')
end
```
### Timelocked trials, rebound=red
![](./images/rebound_timelocked.png)  

### Timelocked trials, rebound=red, baseline=blue  
![](./images/rebound_and_baseline_timelocked.png)

###Fourier analyses
Here we make Fourier decompositions of the time courses in the cropped data for each of the times of interest. Here we use _ft_freqanalysis_ similar to how we earlier calculated PSD and TFR. But rather than focusing on the entire spectral range we confine _ft_freqanalysis_ to only the 16 Hz. This (sort of) represents the primary band where we saw the beta rebound in the TFR plots. Also, we will not get the spectral power of the signal, but keep the complex-valued Fourier representation.

The complex-valued Fourier representation will be used to calculate the cross-spectral density later. Note that we compute three different Fourier decompositions: 
1. One for the beta rebound time of interest.
2. One for the baseline time of interest.
3. One for the beta rebound and baseline combined.

```{r, engine='octave', eval=FALSE}
%% fourier, decomposition Beta rebound

fouriers_rebound  = cell(1, n_events);
fouriers_baseline = cell(1, n_events);
fouriers_combined = cell(1, n_events);
beta_foi = [16 16];

for event_index = 1:n_events
    
    cfg = [];
    cfg.method     = 'mtmfft';
    cfg.output     = 'fourier';
    cfg.taper      = 'hanning';
    cfg.channel    = channels;
    cfg.foilim     = beta_foi;
    cfg.keeptrials = 'yes';
    cfg.pad        = 'nextpow2';
    
    fouriers_rebound{event_index}  = ft_freqanalysis(cfg, tois_rebound{event_index});
    fouriers_baseline{event_index} = ft_freqanalysis(cfg, tois_baseline{event_index});
    fouriers_combined{event_index} = ft_freqanalysis(cfg, tois_combined{event_index});
    
end
```

### Example plot Fourier
Take a look at the data structure. Here we plot the real value of the complex Fourier representation across trials:

```{r, engine='octave', eval=FALSE}
%% example plot Fourier

event_index = 1;
fourier_rebound = fouriers_rebound{event_index};
fourier_baseline = fouriers_baseline{event_index};

cfg = [];
cfg.method = 'svd';

cmb_fourier_rebound  = ft_combineplanar(cfg, fourier_rebound);
cmb_fourier_baseline = ft_combineplanar(cfg, fourier_baseline);

channel_of_interest_tactile = 'MEG0422+0423';
channel_of_interest_occipital = 'MEG2142+2143';
channel_index_tactile = strcmp(cmb_fourier_rebound.label, channel_of_interest_tactile);
channel_index_occipital = strcmp(cmb_fourier_rebound.label, channel_of_interest_occipital);
n_trials = length(cmb_fourier_rebound.trialinfo);

figure('units', 'normalized', 'outerposition', [0 0 1 1]); %% make full screen figure
hold on %% allows for multiple plot calls
% plot all channels
plot(1:n_trials, real(cmb_fourier_rebound.fourierspctrm(:, :)).^2, 'r'); %% plot rebound in red
plot(1:n_trials, real(cmb_fourier_baseline.fourierspctrm(:, :)).^2, 'b'); %% plot baseline in blue
xlabel('Trial no')
ylabel('Power of beta band');
% highlight single trials
plot(1:n_trials, real(cmb_fourier_rebound.fourierspctrm(:, channel_index_tactile)).^2, 'k', 'linewidth', 10); %% plot rebound tactile in black
plot(1:n_trials, real(cmb_fourier_rebound.fourierspctrm(:, channel_index_occipital)).^2, 'm', 'linewidth', 10); %% plot rebound occipital in magenta
```

All channels, red=rebound, blue=baseline

![](./images/fourier_example_all.png)

With highlights, black="tactile", magenta='occipital'

![](./images/fourier_example_all_with_highlights.png)

## Create a lead field and grid
Create a grid around the brain and estimate the lead field for each of the grid points in the brain.

The lead field is an important concept, which may appear confusing at first. Just to make it more confusing, it is also sometimes called the forward model (which was the case for the Minimum-Norm Estimate).

For any given source (a grid point inside the brain) it is calculated how each sensor (magnetometer, gradiometer or electrode) sees (how much T, T/m or V would it pick up) a source with unit strength (1 nAm). One might say that it says: "For a given source, _if_ it is active, how _would_ the different sensors see it"  

```{r, engine='octave', eval=FALSE}
%% leadfield beamformer

cfg = [];
cfg.grad = baseline_data.grad; %% magnetometer and gradiometer specification
cfg.elec = baseline_data.elec; %% electrode specification
cfg.headmodel = headmodel; %% headmodel used
cfg.channel = channels;
cfg.grid.resolution = 1; %% resolution of ...
cfg.grid.unit = 'cm'; %% ... 1 cm
cfg.senstype = sensor_type;

leadfield = ft_prepare_leadfield(cfg);

cfg.channel = 'megmag';
leadfield_mag = ft_prepare_leadfield(cfg);
```

Take a look at the leadfield toghether with the volume conductor model:

```{r, engine='octave', eval=FALSE}
figure('units', 'normalized', 'outerposition', [0 0 1 1]);
hold on
ft_plot_mesh(ft_convert_units(leadfield, 'mm'));
ft_plot_vol(headmodel);
view([-45 20])
```

![](./images/brain_inside_grid.png)

* Q: What does the source model look like? How many sources are there? The plot command below might be of help:

```{r, engine='octave', eval=FALSE}
lf = ft_convert_units(leadfield, 'mm');
lf.pos = lf.pos(lf.inside,:);
ft_plot_mesh(lf,'vertexcolor','r')
```


### Show lead field (code to generate not shown)
1. The center of the blue circle is the grid point of the source  (the size is just for visibility)  
2. The hotter the red color, the more strongly the given sensor would see it  

![](./images/leadfields.png)

Now that we have the Fourier transformed data, the source space, and the lead field, we are finally able to do the source analysis using the DICS beamformer.

## Source analysis using the DICS beamformer
The beamformer is a spatially adaptive filter, allowing us to estimate the amount of activity at any given location in the brain. The inverse filter is based on minimizing the source power (or variance) at a given location, subject to unit-gain constraint. Unit-gain constraint means that, if a source had the power of amplitude 1 and was projected to the sensors by the lead field, the inverse filter applied to the sensors should then reconstruct power of amplitude one at that location. 

In this step, _ft_sourceanalysis_ will take all the previous ingredients and do the source reconstruction. We specify that we want to use DICS with _cfg.method = 'dics'_. In addition, we will give parameters that tell how the DICS beamformer will regularise noise, and ask to keep the spatial filters that is estimated (_cfg.dics.keepfilter = 'yes'_), which we will plot for educational purposes.

We will start by calculating a spatial filter that is shared between the rebound and the baseline based on the combined Fourier analysis. The combined filter is then used on the beta rebound and baseline data separately. 

```{r, engine='octave', eval=FALSE}
%% source analysis

beamformers_rebound  = cell(1, n_events);
beamformers_baseline = cell(1, n_events);
beamformers_combined = cell(1, n_events);

for event_index = 1:n_events
    
    cfg = [];
    cfg.method              = 'dics';                   % Dynamic Imaging of Coherent Sources
    cfg.frequency           = fouriers_rebound{1}.freq; % the frequency from the fourier analysis (as defined above to be 16 Hz)
    cfg.grid                = leadfield;                % Our grid and the leadfield
    cfg.headmodel           = headmodel;                % our headmodel (tells us how the magnetic field/electrical potential is propagated)
    cfg.dics.projectnoise   = 'yes';                    % estimate noise
    cfg.dics.lambda         = '10%';                    % how to regularise
    cfg.dics.keepfilter     = 'yes';                    % keep the spatial filter in the output
    cfg.dics.realfilter     = 'yes';                    % retain the real values
    cfg.channel             = channels;
    cfg.senstype            = sensor_type;
    cfg.grad                = baseline_data.grad;
    cfg.elec                = baseline_data.elec;
    
    beamformers_combined{event_index} = ft_sourceanalysis(cfg, fouriers_combined{event_index});
    
    % Copy filter to cfg
    cfg.grid.filter = beamformers_combined{event_index}.avg.filter;
    
    beamformers_rebound{event_index} = ft_sourceanalysis(cfg, fouriers_rebound{event_index});
    beamformers_baseline{event_index} = ft_sourceanalysis(cfg, fouriers_baseline{event_index});
    
end
```

* Q: Why do you think we use a common filter for both the beta rebound data and baseline data?

### Plot sensitivity profiles of the two channels on the brain

Let us plot the filter. Note the center of the head bias.

```{r, engine='octave', eval=FALSE}
%% plot sensitivity profiles of the two channels on the brain

load mri_resliced.mat

close all

beam = beamformers_rebound{1}; %% choose a given source reconstruction
filter = beam.avg.filter(:); %% get the spatial filter
n_grid_points = length(filter); 
n_channels = length(fouriers_rebound{1}.label);
norms = zeros(n_channels, n_grid_points); %% prepare to get the magnitudes in each of the three directions x, y, z

for n_grid_point = 1:n_grid_points
    for n_channel = 1:n_channels
        if ~isempty(filter{n_grid_point})
            norms(n_channel, n_grid_point) = norm(filter{n_grid_point}(:, n_channel)); %% get the magnitudes for each combination of channels and grid points that are inside
        end
    end
end

channel_of_interest_tactile   = 'MEG0432'; %% a central sensor
channel_of_interest_occipital = 'MEG2142'; %% an occipital sensor
channel_index_tactile   = strcmp(fourier_rebound.label, channel_of_interest_tactile);
channel_index_occipital = strcmp(fourier_rebound.label, channel_of_interest_occipital);

beam.tactile_filter = norms(channel_index_tactile, :)'; %% add field that can be plotted for each channel
beam.occipital_filter = norms(channel_index_occipital, :)';

% interpolate onto MRI

cfg = [];
cfg.downsample = 2;
cfg.parameter = 'tactile_filter';

beam_inter_tactile = ft_sourceinterpolate(cfg, beam, mri_resliced);

cfg = [];
cfg.downsample = 2;
cfg.parameter = 'occipital_filter';

beam_inter_occipital = ft_sourceinterpolate(cfg, beam, mri_resliced);

cfg = [];
cfg.opacitylim = [0 5];
cfg.funparameter = 'tactile_filter';
cfg.location = [-27.5 6.5 100.5];

ft_sourceplot(cfg, beam_inter_tactile);
set(gcf,'units','normalized','outerposition',[0 0 1 1])

cfg = [];
cfg.opacitylim = [0 5];
cfg.funparameter = 'occipital_filter';
cfg.location = [-15.5 -51.5 40.5];

ft_sourceplot(cfg, beam_inter_occipital);
set(gcf,'units','normalized','outerposition',[0 0 1 1])
```

### Tactile channel  
![](./images/tactile_filter.png)

### Occipital channel 
![](./images/occipital_filter.png)

## Plot Source reconstructions

Before plotting the result of the source reconstruction, make an educated guess about where you will see active sources.

```{r, engine='octave', eval=FALSE}
%% plot non-contrasted

cfg = [];
cfg.downsample = 2;
cfg.parameter = 'pow';

beam_inter = ft_sourceinterpolate(cfg, beam, mri_resliced);

cfg = [];
cfg.funparameter = 'pow';

ft_sourceplot(cfg, beam_inter);
set(gcf,'units','normalized','outerposition',[0 0 1 1])
```

![](./images/centre_of_the_head_bias.png)

Do this look like what we should expect? The center of the head bias means that the source reconstructions in themselves are not readily interpretable.  

## Make a contrast of the rebound and the baseline
Here we take the difference between the source reconstructions and divide it by the baseline power. This will give us the relative increase in power in the beta rebound time-window.

```{r, engine='octave', eval=FALSE}
%% contrasts

beamformers_contrasts = cell(1, n_events);

for event_index = 1:n_events
    
    contrast = beamformers_rebound{event_index}; %% just make a copy
    contrast.avg.pow = (beamformers_rebound{event_index}.avg.pow - beamformers_baseline{event_index}.avg.pow) ./ ...
        beamformers_baseline{event_index}.avg.pow;
    
    beamformers_contrasts{event_index} = contrast;
    
end
```

Interpolate pow data onto the MRI of the participant to plot the results:

```{r, engine='octave', eval=FALSE}
%% contrasts interpolated onto MRI

load mri_resliced.mat
beamformers_contrasts_interpolated = cell(1, n_events);

for event_index = 1:n_events
    
    cfg = [];
    cfg.downsample = 2;
    cfg.parameter = 'pow';
    
    beamformers_contrasts_interpolated{event_index} = ft_sourceinterpolate(cfg, beamformers_contrasts{event_index}, mri_resliced);
    
end
```

Now plot all the contrasts:

```{r, engine='octave', eval=FALSE}
%% plot source analyses

close all
names = {'right_little_finger' 'right_ring_finger' 'right_middle_finger' 'right_index_finger' 'right_thumb'};

for event_index = 1:n_events
    
    to_plot = beamformers_contrasts_interpolated{event_index};
    max_pow = max(to_plot.pow);
    min_pow = min(to_plot.pow);

    cfg = [];
    cfg.method = 'ortho';
    cfg.funparameter = 'pow';
    cfg.maskparameter = cfg.funparameter;
    cfg.funcolorlim = [0 max_pow];
    cfg.opacitylim = [0 max_pow];
    cfg.location = [-29.5 10.5 100.5];

    ft_sourceplot(cfg, to_plot);
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    
end
```

### Right little finger
![](./images/beamformer_contrast_right_little_finger.png)

### Right ring finger
![](./images/beamformer_contrast_right_ring_finger.png)

### Right middle finger
![](./images/beamformer_contrast_right_middle_finger.png)

### Right index finger
![](./images/beamformer_contrast_right_index_finger.png)

### Right thumb
![](./images/beamformer_contrast_right_thumb.png)


* Q: The source reconstruction does not look very good for the index finger. Give an explanation why this might be the case (Hint: Look at the TFR responses).
