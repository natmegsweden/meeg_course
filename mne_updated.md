---
title: "Minimum-Norm Estimate (MNE)"
author: "Mikkel Vinding & Lau Møller Andersen"
date: "Feb. 2018; NatMEG, Stockholm, Sweden"
output: html_document
---
  
## Set up paths
Change these to appropriate paths for your operating system and setup

```{r, engine='octave', eval=FALSE}
addpath /home/lau/matlab/fieldtrip-20170613_workshop_natmeg/
ft_defaults

raw_meg_path = '/archive/20067_workshop_source_reconstruction/MEG/';
meg_path = '/home/share/workshop_source_reconstruction/data/MEG/';
mri_path = '/home/share/workshop_source_reconstruction/data/MRI/';
```

### Set up subject-specific paths
Make subject and recording specific paths (the cell array "subjects_and_dates" can be expanded)

```{r, engine='octave', eval=FALSE}
%% subjects and dates

subjects_and_dates = ...
{
  'NatMEG_0177/170424/'
};

output_path = fullfile(meg_path, subjects_and_dates{1});
cd(output_path)

events = [1 2 4 8 16]; %% right little, right ring, right middle, right index, right thum

```

## Import source space
Minimum-norm estimate (MNE) use a cortical surface as source model. We will not run the preparation of the cortical surface as it takes about 10 hours. You can find tutorial documentation on how the cortical surface source model was made here: http://natmeg.se/onewebmedia/prepare_MNE_sourcespace.html. In the tutorial, you can read how to export MRI to Freesurfer to extract the cortical surface, and set up a source model with MNE-C. This creates a set of point equally sampled across the cortical surface. These points are our source model for MNE. 

You will find the file "Sub02-oct-6-src.fif" in the tutorial data files.

We read the MNE-C source space into MATLAB by using _ft_read_headshape_ from the MRI folder:
  
```{r, engine='octave', eval=FALSE}
cd(fullfile(mri_path,'freesurfer/Sub02/bem'))  % This should change if looping over subjects

sourcemodel = ft_read_headshape('Sub02-oct-6-src.fif', 'format', 'mne_source'); 

cd(fullfile(meg_path,subjects_and_dates{1}))
````

* Q: What does the "sourcemodel" structure contain?

Note that the units are in meters. You can use ft_convert_units to make the units milimeters instead:
  
```{r, engine='octave', eval=FALSE}
cd(fullfile(/home/share/workshop_source_reconstruction/data/MRI/))
sourcemodel = ft_convert_units(sourcemodel, 'mm');
````

You can plot the source model with the function ft_plot_mesh:
  
```{r, engine='octave', eval=FALSE}
figure; ft_plot_mesh(sourcemodel, 'edgecolor', 'k'); camlight 
````

![](./images/mne_sourcespace.png)

Looks nice!
  
Each point on the grid will represent a dipole in the MNE source reconstruction.

* Q: How many dipoles will the MNE source reconstruction contain?
  
### Load data
```{r, engine='octave', eval=FALSE}
load headmodel_mne_meg.mat
load timelockeds.mat

````

Now you should have the three main ingredients for doing MNE source reconstruction with MEG:
  
1. Evoked MEG data.
2. A volume model of the brain.  
3. A source space model of the cortical surface.  

Take a look at the headmodel and sourcemodel toghether:
  
```{r, engine='octave', eval=FALSE}
figure; hold on
ft_plot_vol(headmodel_mne_meg, 'facealpha', 0.5, 'edgecolor', 'none');
ft_plot_mesh(sourcemodel, 'edgecolor', 'k'); camlight 

````

![](./images/mne_sourcespace2.png)

* Q: What is going on -- explain the picture?
  
### Aligning source space and volume model
  
We have aligned the MRI that we used to create the volume model to the Neuromag coordinate system used by the MEG scanner. However, to make the cortical surface in Freesurfer, we had to convert the MRI to the spm coordinate system and exported it to Freesurfer (see http://natmeg.se/onewebmedia/prepare_MNE_sourcespace.html). When we imported the source space back into MATLAB, it kept this coordinate system. Our head model and source space are now in two different coordinate systems.

To get the cortical surface back into the head, you need to transform the points in source space to the coordinate system of the head model. Load the transformation matrix of the resliced MRI in the spm coordinate system (_transform_vox2spm_rs_). We will use this to create a transformation from the spm coordinate system to the neuromag coordinate system. If you have run the prepare MNE source space tutorial, you must read the transformation files you created yourself.


```{r, engine='octave', eval=FALSE}
% from voxel to spm-coordsys.
load transform_vox2spm.mat

% Load segmented mri and get transform
load mri_segmented_mne.mat

% Get transformation matrix from voxel to neuromag-coordsys
transform_vox2neuromag = mri_segmented_mne.transform;

% Get transformation matrix from spm to neuromag
T = transform_vox2neuromag/transform_vox2spm_rs;
````

Now transform the source space using the transformation matrix **T**. Note, that if you run the following line several times, it will apply the transformation each time, transforming from the current position to a new location, which will be wrong. If in doubt, save the source model before proceeding.

```{r, engine='octave', eval=FALSE}
sourcemodelT = ft_transform_geometry(T, sourcemodel);
````

Be aware that each time you run _ft_transform_geometry_ it will apply the transform **T** to the coordinated of the points in the source model.

Take a look at the headmodel and sourcespace again:
  
```{r, engine='octave', eval=FALSE}
figure; hold on
ft_plot_vol(headmodel_mne_meg, 'facealpha', 0.5, 'edgecolor', 'none');
ft_plot_mesh(sourcemodelT, 'edgecolor', 'k'); camlight 

````

![](./images/mne_sourcespace3.png)

Remember to save the source model.

```{r, engine='octave', eval=FALSE}
save('sourcemodelT.mat','sourcemodelT')
disp('Done');
````

## MNE source reconstruction on MEG data

Now we can create the leadfield and do the source reconstruction. Load the relevant files (if you do not already have then loaded):
  
```{r, engine='octave', eval=FALSE}
load timelockeds.mat
load sourcemodelT.mat
load headmodel_mne_meg.mat
disp('done')
````

For this tutorial we will use the stimulation of the indexfinger -- corresponding the the 4th condition:
  
```{r, engine='octave', eval=FALSE}
data_meg = timelockeds{4};
````

Calculate the leadfield using ft_prepare_leadfield with the MEG data, the source model, and the appropriate head model. Here we use the gradiometers (meggrad).

```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.grad                = data_meg.grad;              % sensor positions
cfg.channel             = 'meggrad';                  % the used channels
cfg.senstype            = 'meg';
cfg.grid.pos            = sourcemodelT.pos;           % source points
cfg.grid.inside         = 1:size(sourcemodelT.pos,1); % all source points are inside of the brain
cfg.headmodel           = headmodel_mne_meg;              % volume conduction model

leadfield_mne = ft_prepare_leadfield(cfg,data_meg);
````

Take a look at the "leadfield_mne" structure.

* Q (optional): How does it compare to the leadfield structure you created for the dipole tutorial (hint: look at the size of the "pos" matrix or try to plot the points)?

Plot the point in the source model. Require that you have the leadfield from the dipole fit tutorial loaded in memory. Otherwise you have to re-load it (or skip this part).

```{r, engine='octave', eval=FALSE}
plot3(leadfield_mne.pos(:,1),leadfield_mne.pos(:,2),leadfield_mne.pos(:,3),'o')  % o's not zero's
plot3(leadfield.pos(:,1),leadfield.pos(:,2),leadfield.pos(:,3),'o') % o's not zero's
````

## MNE source reconstruction
Finally, we do the source reconstruction using _ft_sourceanalysis_. Specify _cfg.method = 'mne'_ to use MNE.

MNE require additional parameters (lambda) that indicate how to scale the noise covariance. Note that the noise covariance matrix already is in the "data" structure. This was calculated already in the preprocessing steps when we calculated the evoked fields. 

```{r, engine='octave', eval=FALSE}
cfg                     = [];
cfg.method              = 'mne';
cfg.channel             = 'meggrad';
cfg.senstype            = 'meg';
cfg.grid                = leadfield_mne;
cfg.headmodel           = headmodel_mne_meg;
cfg.mne.prewhiten       = 'yes';
cfg.mne.lambda          = 3;
cfg.mne.scalesourcecov  = 'yes';
source_mne  = ft_sourceanalysis(cfg,data_meg);

save source_mne.mat source_mne;
disp('Done')
````

Look at the _source_mne_ structure -- especially the _source_mne.avg_ structure. 

* Q: What is in the output structure -- what is the dimension of the data? 
  
## Visualize
  
Plot the MNE source reconstruction on the grid. The source_mne structure contains values for all grid points and all time points. To visualize the source reconstruction on the grid, we need to decide what time points to plot.

Go back to the evoked data to find some interesting times to plot using ft_multiplotER. Since we used the gradiometers, we will start by plotting the fields from the combined gradiometers. For example, let us compare the first peak around 50 ms and the large peak from 140-150 ms. Use the drag and click tool on the multiplot to plot topographical plots.

```{r, engine='octave', eval=FALSE}
data_cmb = ft_combineplanar([],data_meg); %Combine gradiometers

cfg = [];
cfg.layout = 'neuromag306cmb.lay';
figure; ft_multiplotER(cfg, data_cmb);
````

Find peaks in sensor space, look at the topographies, and try to guess what the sources might look like?
  
![](./images/peak_mne1.png)
![](./images/peak_mne2.png)
![](./images/peak_mne3.png)
![](./images/peak_mne4.png)
![](./images/mne_gradTopo.png)

Select the time-windows of interest and plot on the sourcespace:
  
```{r, engine='octave', eval=FALSE}
source_mne.tri = sourcemodelT.tri; %ft_sourceplot need the triangulation information, add to source reconstruction.

cfg = [];
cfg.method          = 'surface';
cfg.funparameter    = 'pow';
cfg.funcolormap     = 'jet';
cfg.latency         = .060;     % The time-point to plot
cfg.colorbar        = 'no';

ft_sourceplot(cfg, source_mne)
````

![](./images/mne_50ms.png)

* Q. By a glance, where do you see activated sources?  
* Q: How many different cortical patches are "active" at this time point?
  
This image shows a single time point for all estimated sources. Each source (i.e., each grid-point on the cortical surface) has a time series of activation. You can in principle treat each source as its own "channel" -- you can view the activation of all sources over time, analogous to the activation in MEG or EEG channels. Use MATLAB to plot source activation across time (here only plotting every third time-series to save memory)

```{r, engine='octave', eval=FALSE}
plot(source_mne.time,source_mne.avg.pow(1:3:end,:));
````

![](./images/mne_all_sources.png)

Now try to use _ft_sourceplot_ to plot the source space topography of other time points.

* Q: Can you find the time-point corresponding to the image below? Try to look at the scalp topographies to get an educated guess about the time.  

![](./images/mne_150ms.png)

Try to visualize the sources of some of the following components from the scalp data. You can visualize the average over time over with _ft_sourceplot_ by changing the _cfg.latency_ to an interval [start stop], and give an additional cfg argument _cfg.avgovertime = 'yes'_.


```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.method          = 'surface';
cfg.funparameter    = 'pow';
cfg.funcolormap     = 'jet';
cfg.latency         = [.350 .40]; % average actvity from 350 ms to 400 ms
cfg.avgovertime     = 'yes';
cfg.colorbar        = 'no';

ft_sourceplot(cfg, source_mne)
````

![](./images/mne_350ms.png)


To get an overview of both the temporal and spatial features of the source reconstruction, you can use _ft_sourcemovie_ to make a movie of the source activation. Note that the triangulation field must be present in the _source_mne_ structure for _ft_sourcemovie_ to work. You can play around with the scaling to get a nicer looking movie.

```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.funparameter = 'pow';
ft_sourcemovie(cfg,source_mne);
````

![](./images/mne_movie.png)

It is difficult to answer what type of source model and source reconstruction method that is better. The matter model depends on the question you want to answer and the kind of signal you are interested in for the particular answer. For example, if you know that a focal patch of cortex generates the signal then a single dipole model might be sufficient, and you do not gain extra information by including the entire cortex. Similar, if we know that the process we are interested in requires distributed sources then itis not valid to assume that a single dipole is sufficient.

* Q: Around what times in the time-series might a dipole model be sufficient, and at what points is a distributed model better?