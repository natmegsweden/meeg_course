---
title: "Dipole Fitting"
author: "Lau Møller Andersen & Mikkel Vinding"
date: "16 June 2017; NatMEG, Stockholm, Sweden"
output: html_document
---

## Please read this about paths and files before you proceed!

Note that if work from the path where all the downloadable parts are downloaded to, you don't need to change the paths  
(leave them as is, but do evaluate the sections)  
The paths simply serve as an example of how you can set up your analysis structure, which is especially useful if you have more than one subject 

### Set up general paths
Change these to appropriate paths for your operating system and setup  


```{r, engine='octave', eval=FALSE}
%% paths

addpath /home/lau/matlab/fieldtrip-20170613_workshop_natmeg/
ft_defaults

raw_meg_path = '/archive/20067_workshop_source_reconstruction/MEG/';
meg_path = '/home/share/workshop_source_reconstruction/data/MEG/';
mri_path = '/home/share/workshop_source_reconstruction/data/MRI/';
set(0, 'DefaultAxesFontWeight', 'bold', 'defaultaxesfontsize', 20, 'defaultlinelinewidth', 2);
```

### Set up subject specific paths
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

### Load relevant files

We are loading three different files here:  

1. We are loading the timelocked data for the MEG and EEG data
2. We are loading the headmodels for the MEG data
3. We are loading the headmodels for the EEG data
  
(For now leave the tissue_based_headmodel on as 'yes', later we will play around with it)

```{r, engine='octave', eval=FALSE}
%% go to relevant path and load data

tissue_based_headmodel = 'yes';

% cd(output_path) %% comment this in if you want to set the output path
disp 'Loading timelockeds and headmodel'
load timelockeds.mat
load headmodel_meg.mat
load headmodel_eeg.mat
load headmodel_eeg_fourspheres_elec.mat
load headmodel_meg_singlesphere_headshape.mat
disp Done

if ~strcmp(tissue_based_headmodel, 'yes');
    headmodel_eeg = headmodel_eeg_fourspheres_elec;
    headmodel_eeg.type = 'concentricspheres';
    headmodel_meg = headmodel_meg_singlesphere_headshape;
end
```

### Identify components of interest

```{r, engine='octave', eval=FALSE}
%% identify components of interest

close all
figure('units', 'normalized', 'outerposition', [0 0 1 1]);

% multiplot over all magnetometers
cfg = [];
cfg.layout = 'neuromag306mag.lay'; % layour for magnetometers

ft_multiplotER(cfg, timelockeds{4}); %% explore this by yourself

% combine planar gradiometers (makes them interpretable)

cfg = [];

cmb = ft_combineplanar(cfg, timelockeds{4});

% plot channels "over" SI

cfg = [];
cfg.channel = {'MEG0412+0413' 'MEG0422+0423' 'MEG0432+0433' 'MEG0442+0443'}; %% channels "over" SI
cfg.layout = 'neuromag306cmb.lay'; %% layout for combined gradiometers

figure('units', 'normalized', 'outerposition', [0 0 1 1]);

ft_singleplotER(cfg, cmb);
xlabel('Time (s)')
ylabel('Root mean squared activity')

% plot topoplot for early activity
cfg = [];
cfg.layout = 'neuromag306cmb.lay';
cfg.xlim = [0.045 0.065]; % s
cfg.comment = 'no';

figure('units', 'normalized', 'outerposition', [0 0 1 1]);

ft_topoplotER(cfg, cmb);

% plot topoplot for later activity

cfg.xlim = [0.115 0.155]; % s
cfg.zlim = [0 6e-12]; % scale so both peaks can be seen

figure('units', 'normalized', 'outerposition', [0 0 1 1]);

ft_topoplotER(cfg, cmb);
``` 
### The evoked time courses for the "SI" sensors  
Notice the peaks around 55 msec and 135 msec

![](./images/evoked_cmb_SI_and_SII.png)

### A topographical plot averaging over 45-65 msec over SI


![](./images/topp_cmb_SI.png)

### A topographical plot averaging over 115-155 msec over SII

![](./images/topo_cmb_SII.png)

### Create a leadfield and grid
Create a grid around the brain and estimate the leadfield for each of the grid points in the brain  

#### The leadfield is an important concept, which may appear confusing at first  
1. For any given source (a grid point inside the brain) it is calculated how each sensor (magnetomer, gradiometer or electrode) sees (how much T, T/m or V would it pick up) a source with unit strength (1 nAm)  
2. One might say that it says: "For a given source, _if_ it is active, how _would_ it be seen by the different sensors"  
3. It is also sometimes called the forward model  

We define the leadfields for each of the three kinds of data here (magnetometers, gradiometers and electrodes)

```{r, engine='octave', eval=FALSE}
%% make leadfields eeg and meg

cfg = [];
cfg.headmodel = headmodel_eeg;
cfg.elec = timelockeds{1}.elec;
cfg.senstype = 'eeg';
cfg.grid.resolution = 1;
cfg.grid.unit = 'cm';

leadfield_eeg = ft_prepare_leadfield(cfg, timelockeds{1});

cfg.senstype = 'meg';
cfg.grad = timelockeds{1}.grad;
cfg.headmodel = headmodel_meg;
cfg.channel = 'megmag';

leadfield_mag = ft_prepare_leadfield(cfg, timelockeds{1});

cfg.channel = 'meggrad';

leadfield_grad = ft_prepare_leadfield(cfg, timelockeds{1});
```

```{r, engine='octave', eval=FALSE}
%% plot grid and headmodel

figure('units', 'normalized', 'outerposition', [0 0 1 1]);
hold on
ft_plot_mesh(ft_convert_units(leadfield_mag, 'mm'));
ft_plot_vol(headmodel_meg);
view([-45 20]) %% try to rotate it yourself with the plot tools
```

### Show brain inside grid
![](./images/brain_inside_grid.png)

### Show leadfield (code to generate not shown)
1. The centre of the blue circle is the grid point of the source  (the size is just for visibility)  
2. The hotter the red colour, the more strongly the given sensor would see it  

![](./images/leadfields.png)

### Fit dipoles

1. We fit single dipoles separately for the three kinds of sensors (magnetometers, gradiometers and electrodes)  
2. We do it for two latencies.  
  2a. Early: 45-65 msec, probably contralateral SI  
  2b. Late: 115-155 msec, probably bilateral SII  
3. Our optimization parameter is residual variance which we try to minimize. This is how much data that is left unexplained by the dipole model  

```{r, engine='octave', eval=FALSE}
%% dipole fits

n_events = length(events); % the events to be looped through

% the six cell arrays below are preparation for the fits

dipoles_mag_early =  cell(1, n_events);
dipoles_grad_early = cell(1, n_events);
dipoles_eeg_early =  cell(1, n_events);

dipoles_mag_late =  cell(1, n_events);
dipoles_grad_late = cell(1, n_events);
dipoles_eeg_late =  cell(1, n_events);

early_latency = [0.045 0.065]; % s
late_latency = [0.115 0.155]; % s

for event_index = 1:n_events

    cfg = [];
    cfg.gridsearch = 'yes'; %% search the grid for an optimal starting point
    cfg.grid = leadfield_mag; %% supply the grid/leadfield
    cfg.headmodel = headmodel_meg; %% supply the headmodel
    cfg.dipfit.metric = 'rv'; %% the metric to minimize (the relative residual variance: proportion of variance left unexplained by the dipole model)
    cfg.model = 'regional'; %% we assume that the dipole has a fixed position during the time points in the latency range
    cfg.senstype = 'meg'; %% sensor type
    cfg.channel = 'megmag'; %% which channels to use
    cfg.nonlinear = 'yes'; %% do a non-linear search

    % magnetometer fits
    cfg.latency = early_latency; %% specify the latency
    cfg.numdipoles = 1; %% we only expect contralateral activity
    cfg.symmetry = []; %% empty for single dipole fit
    dipoles_mag_early{event_index} = ft_dipolefitting(cfg, timelockeds{event_index});
    
    cfg.latency = late_latency;
    cfg.numdipoles = 2; %% we expect bilateral activity
    cfg.symmetry = 'x'; %% we expect it to be symmetrical in the x-direction (ear-to-ear)
    dipoles_mag_late{event_index} = ft_dipolefitting(cfg, timelockeds{event_index});
 
    % gradiometer fits
    cfg.channel = 'meggrad';
    cfg.grid = leadfield_grad;

    cfg.latency = early_latency;
    cfg.numdipoles = 1;
    cfg.symmetry = [];
    dipoles_grad_early{event_index} = ft_dipolefitting(cfg, timelockeds{event_index});
    
    cfg.latency = late_latency;
    cfg.numdipoles = 2;
    cfg.symmetry = 'x';
    dipoles_grad_late{event_index} = ft_dipolefitting(cfg, timelockeds{event_index});
    
    %% eeg fits
    cfg.senstype = 'eeg';
    cfg.channel = 'eeg';
    cfg.grid = leadfield_eeg;
    cfg.headmodel = headmodel_eeg;
    
    cfg.latency = early_latency;
    cfg.numdipoles = 1;
    cfg.symmetry = [];
    dipoles_eeg_early{event_index} = ft_dipolefitting(cfg, timelockeds{event_index});
    
    cfg.latency = late_latency;
    cfg.numdipoles = 2;
    cfg.symmetry = 'x';
    dipoles_eeg_late{event_index} = ft_dipolefitting(cfg, timelockeds{event_index});
    
end

disp('Done')
```

### Plot EARLY dipoles

Try to rotate these with the normal plot tools to get a feeling of where the dipole is located  

```{r, engine='octave', eval=FALSE}
%% plot dipoles_early on brain

close all
load mri_resliced_cm.mat

all_dipoles = {dipoles_mag_early dipoles_grad_early dipoles_eeg_early};
colours = {'r' 'g' 'y' 'b' 'm'};
fingers = {'right little finger' 'right ring finger' 'right middle finger' 'right index finger' 'right thumb'};
% figure('units', 'normalized', 'outerposition', [0 0 1 1]);

for dipole_type_index = 1:length(all_dipoles)
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    for event_index = 1:n_events
        subplot(2, 3, event_index)
        dipole = all_dipoles{dipole_type_index}{event_index};

        hold on
        ft_plot_dipole(dipole.dip.pos(1, :), mean(dipole.dip.mom(1:3, :), 2), 'color', colours{event_index});

        pos = mean(dipole.dip.pos, 1);

        ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1);
        ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1);
        ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1);

        title_text = sprintf([fingers{event_index} ' at coordinates (cm):\nx: ' num2str(pos(1)) '; y: ' num2str(pos(2)) '; z: ' num2str(pos(3))]);
        title(title_text, 'fontsize', 24)

        axis tight
        axis off

    end

end
```

### Plot LATE dipoles

Try to rotate these with the normal plot tools to get a feeling of where the dipoles are located  

```{r, engine='octave', eval=FALSE}
%% plot dipoles_late on brain

close all
load mri_resliced_cm.mat
% load dipoles_late.mat

all_dipoles = {dipoles_mag_late dipoles_grad_late dipoles_eeg_late};
colours = {'r' 'g' 'y' 'b' 'm'};
fingers = {'right little finger' 'right ring finger' 'right middle finger' 'right index finger' 'right thumb'};

for dipole_type_index = 1:length(all_dipoles)
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    for event_index = 1:n_events
        subplot(2, 3, event_index)
        dipole = all_dipoles{dipole_type_index}{event_index};

        hold on
        ft_plot_dipole(dipole.dip.pos(1, :), mean(dipole.dip.mom(1:3, :), 2), 'color', colours{event_index});
        ft_plot_dipole(dipole.dip.pos(2, :), mean(dipole.dip.mom(4:6, :), 2), 'color', colours{event_index});

        pos = mean(dipole.dip.pos, 1);

        ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1);
        ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1);
        ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1);

        title_text = sprintf([fingers{event_index} ' at coordinates (cm):\nx: ' num2str(pos(1)) '; y: ' num2str(pos(2)) '; z: ' num2str(pos(3))]);
        title(title_text, 'fontsize', 24)

        axis tight
        axis off

    end
end
```

### Dipole EARLY plot magnetometers
We get a good fit to the SI

![](./images/dipoles_early_2magnetometers.png)

### Dipole EARLY plot gradiometers
We get a good fit to the SI

![](./images/dipoles_early_2gradiometers.png)

### Dipole EARLY plot electrodes
We get a bad fit to the SI

![](./images/dipoles_early_2electrodes.png)

### Dipole LATE plot magnetometers
We get a bad fit to the SII

![](./images/dipoles_late_2magnetometers.png)

### Dipole LATE plot gradiometers
We get a good fit to the SII

![](./images/dipoles_late_2gradiometers.png)

### Dipole LATE plot electrodes
We get a bad fit to the SII

![](./images/dipoles_late_2electrodes.png)

## Intermittent conclusions
We get the best fits from the gradiometers. In contrast the electrodes seem way off.  

* Now try and change _tissue_based_headmodel_ to 'no' in the _go to relevant path and load data_ block (the third one). Then rerun the whole script by pressing F5. This will base the analysis on sphere models instead of tissue based models  

### The concentric spheres head model

1. It consists of four spheres  
    1a. Brain  
    1b. Cerebro-spinal fluid  
    1c. Skull  
    1d. Skin  
2. They have the following radii and conductivities taken from the BESA software:
     radius       [71 72 79 85] cm
     conductivity [0.3300 1 0.0042 0.3300]
3. It will be used for the EEG data

![](./images/concentric_head_model.png)

### The single sphere headmodel

It consists of a single sphere, which model the brain (the magnetic field does not have different conductivities through the different media)  

![](./images/single_sphere_head_model.png)

### Dipole EARLY plot magnetometers _single sphere_
We get a good fit to the SI

![](./images/dipoles_early_2_concentricspheres_magnetometers.png)

### Dipole EARLY plot gradiometers _single sphere_
We get a good fit to the SI

![](./images/dipoles_early_2_concentricspheres_gradiometers.png)

### Dipole EARLY plot electrodes _concentric spheres_
_Now_ We get a good bad fit to the SI

![](./images/dipoles_early_2_concentricspheres_electrodes.png)

### Dipole LATE plot magnetometers _single sphere_
We get a bad fit to the SII

![](./images/dipoles_late_2_concentricspheres_magnetometers.png)

### Dipole LATE plot gradiometers _single sphere_
We get a good fit to the SII

![](./images/dipoles_late_2_concentricspheres_gradiometers.png)

### Dipole LATE plot electrodes _concentric spheres_
_Now_ we get a good fit to the SII

![](./images/dipoles_late_2_concentricspheres_electrodes.png)

## Conclusions

1. With the non-tissue related models, we can get reasonable fits with the EEG as well, but spherical models do not make sense for the more advanced modelling techniques such as beamformer and minumum norm estimates that we will turn to next.  
2. This underlines that the tissue-based EEG headmodels are hard to create (It may have to do with a poor MRI), but in general the advantage of MEG is that you just need to be able to delineate the brain.
