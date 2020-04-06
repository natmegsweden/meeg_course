---
title: "Preprocess MRI data"
author: "Lau Møller Andersen & Mikkel C. Vinding"
date: "Feb. 2018; NatMEG, Stockholm, Sweden"
output: html_document
---

To do source reconstruction of MEG and EEG signals we need to solve the inverse problem -- i.e., find the sources that generate the magnetic or electric field patterns that we measure. This inverse problem has infinite solutions. To be able to estimate the sources of the magnetic and electrical signals we need constraints on the possible solutions. Lucky enough, we know that the origin of MEG and EEG signals is not any random electric currents, but currents in the brain -- more precise the pyramidal cells in neocortex. We can use this information to constrain the solutions to the inverse problem to a set of pre-specified locations. We assume that the activity we measure comes from the brain and thus limit our possible sources to be within the brain.

We do this by constraining the inverse solution based on the anatomy of the brain and head, information about the volume containing the electric fields (the volume conduction model), and information about the location of the head relative to the MEG/EEG sensors.

To do source reconstruction of MEG and EEG signals, we need a model of the head and how it conducts electrical currents. We obtain a model of the head and the brain from a structural magnetic resonance image (MRI). There three basic ingredients in MEG and EEG source reconstruction is:

1. MEG/EEG data.
2. A structural MRI.
3. Information about the relative location of sensors to the head.

Usually, we will want individual MRI for each participant to accommodate individual differences in head geometry and shape of the brain, but in some cases, it can be sufficient to use a template brain.

In this tutorial, we have an MRI for the subject. The raw dicom files are located in the folder "dicoms". The MEG/EEG data is in the raw .fif files. The information about the relative location of sensors to the head is also contained in the raw fif file. This is the head point we "drew" on the subjects head together with the locations of the EEG electrodes and information about the MEG sensors location relative to the HPI coils, we measured inside the MEG.

In this tutorial, we will create a volume conductor model of the subjects head (one for MEGF and one for EEG), and make sure it is aligned with the position.


### Set up general paths
Change these to appropriate paths for your operating system and setup

```{r, engine='octave', eval=FALSE}
restoredefaultpath
addpath /home/lau/matlab/fieldtrip-20170613_workshop_natmeg/
ft_defaults

raw_meg_path = '/archive/20067_workshop_source_reconstruction/MEG/';
meg_path = '/home/share/workshop_source_reconstruction/data/MEG/';
mri_path = '/home/share/workshop_source_reconstruction/data/MRI/';
```

### Set up subject specific paths
Make subject and recording specific paths (the cell array "subjects_and_dates" can be expanded)

```{r, engine='octave', eval=FALSE}
%% subjects_and_dates

subjects_and_dates = ...
                    {
                        'NatMEG_0177/170424/'
                    };
                
subjects = ...
            {
                'NatMEG_0177'
            };
                
filenames = {...
        'tactile_stim_raw_tsss_mc.fif'
        'tactile_stim_raw_tsss_mc-1.fif'
        'tactile_stim_raw_tsss_mc-2.fif'
            };
        
output_path = fullfile(meg_path, subjects_and_dates{1});
```

## Read in dicom files
Dicom is the standard file format for MRI. We will read the dicom files with the function _ft_read_mri_ to read the MRI volume. To read the dicoms into a single volume, it is enough to provide the filename of any of the dicom in the sequence; _ft_read_mri_ will recognize the rest of the files and read all of them (assuming they are in the same directory).

```{r, engine='octave', eval=FALSE}
dicom_name = '00000001.dcm';

dicom_path = fullfile(mri_path, 'dicoms', subjects{1}, dicom_name);

mri = ft_read_mri(dicom_path);

save([output_path 'mri'], 'mri');
```

Take a look at the newly create _mri_ structure.

You can plot the MRI volume with the function _ft_sourceplot _. Be patient as it might take some time to generate the plots due to the large size of the data. (if it takes too long to plot the volumes, just skip the plotting that is only for inspection purposes).

```{r, engine='octave', eval=FALSE}
figure
ft_sourceplot([], mri);
```

Do not mind that the function is called _ft_sourceplot_ -- it is not yet the source space we are looking at.


![](./images/raw_mri.jpg)

The subject appears to be upside down.

## Co-register MEG/EEG and MRI data
To align the MRI image to the head points in the MEG/EEG data, the first step is to make sure they are in the same coordinate system. As a default, dicom files are in the native coordinate system of the MRI scanner, whereas the MEG and EEG data is in a coordinate system defined by the MEG scanner. These coordinate systems are defined in different ways. 

There a different ways to define head coordinate systems. Different devises use different default coordinate systems. You can read more about how the coordinate systems are defined here: http://www.fieldtriptoolbox.org/faq/how_are_the_different_head_and_mri_coordinate_systems_defined

The tutorial dataset was collected on an Electra Neuromag MEG scanner and use the Neuromag coordinate system. The Neuromag coordinate system is defined form the external landmarks -- i.e. the fiducials, we defined as the first step when we digitized the subject. The coordinate system is defined like this:

* X-axis from the origin towards the RPA point (Right is "up" on X-axis)
* Y-axis from the origin towards the nasion (posterior is "up" on Y-axis)
* Z-axis from the origin upwards orthogonal to the XY-plane (superior is "up" on Z-axis)
* Origin: The intersection of the X-axis and Y-axis.

Note, that through we transform everything to the Neuromag coordinate system, this procedure is the same for any head coordinate system your MEG or EEG data is defined in, whether given by other MEG/EEG device manufacturers. The difference is only what system we use -- usually the "native" system of the MEG/EEG data. Also note, that this is just a practical formality to co-register the MER and MEG/EEG data. Tranfomration to anatomical coordinates such as Talairach-Tournoux or MNI is a different transformation, usually done at a later point (or some cases earlier) of the processing.

### Indicate coordinate system
The first step is to define the axes of the coordinate system. We do this with the FieldTrip function _ft_determine_coordsys_. When you call _dt_determine_coordsys_ a figure will pop up with the axes ontop the MRI. 

```{r, engine='octave', eval=FALSE}
%% determine coordinate system

mri_coordsys = ft_determine_coordsys(mri);

```

![](./images/mri_coordsys.jpg)

In the MATLAB terminal, FieldTrip is asking you to provide information about the direction of the axes. Look at the figure and provide the correct answer to Field Trip. First (r)ight or (l)eft for the X-axis, then (a)nterior or (p)osterior for the Y-axis, and (s)uperior or (?) [!!!] for the Z-axis.

* Q: What is the direction of the axes?

When you have defined the axes, it is a good idea to save the MRI.

```{r, engine='octave', eval=FALSE}
save([output_path 'mri_coordsys'], 'mri_coordsys');
```

### Co-register MRI to Neuromag based on fiducials (nasion and left and right pre-auricular points)
Now we will co-register the MRI image to the head points. In FieldTrip this is done in a step-wise procedure.

First, we align the MRI to the fiducial points. To align MRI and fiducals, we need to specify the fiducials on the MRI. We do this with the FieldTrip function _ft_volumerealign_. This will again plot the volume. You now have to mark the approximate location of the fiducial points NAS, LPA, and RPA.


```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.method   = 'interactive';
cfg.coordsys = 'neuromag';

mri_realigned_1 = ft_volumerealign(cfg, mri_coordsys);

```

Click on the image until the crosshair is on a location corresponding to a fiducial point. You can also scroll through the volume with the arrow keys. Once the crosshair is in a location you think is correct, press "n" for NAS, "l" for LPA, and "r" for RPA on your keyboard. The coordinates of the fiducials will be written in the MATLAB terminal. Press "f" on your keyboard to toggle fiducial visibility on the plot.

You can mark the fiducials multiple times until you are satisfied with the positions. Do it as precise as you can, but do not get upset if it is not perfect. Imprecision will be adjusted in the following steps.

* Once you are satisfied with the positions, save the figure (or take a screenshot) of the location of the fiducials. Make sure the fiducials are visible. Include this is your report.

Finally, you are asked to provide an extra control point that should have a positive z-value point. Click somewhere on the volume that you are sure has a positive Z-value (e.g., the top of the scalp) and press z.

At this point, you might want to save the data again.

### Co-register with the extra digitization points (along the scalp and the face)
Nwo we will load the headpoint from the MEG/EEG data with _ft_read_headshape_. This read the head points from the file. 
```{r, engine='octave', eval=FALSE}
input_path  = fullfile(meg_path, subjects_and_dates{1}, filenames{1});

headshape = ft_read_headshape(input_path);

% Plot head points

[!!!!]
```

![](./images/headpoints.jpg)

We then use _ft_volumerealight_ again, but this time we specify that we want to use an automatic procedure called "iterative closest point" (icp) to fit fiducials and head points.

```{r, engine='octave', eval=FALSE}

cfg = [];
cfg.method              = 'headshape';
cfg.headshape.headshape = headshape;
cfg.headshape.icp       = 'yes';
cfg.coordsys            = 'neuromag';

mri_realigned_2 = ft_volumerealign(cfg, mri_realigned_1);

```

![](./images/headalign.jpg)

### Check co-registration
This step is really just used to check whether the coregistration went okay.
```{r, engine='octave', eval=FALSE}
%% checking co-registration
cfg = [];
cfg.method = 'headshape';
cfg.headshape.headshape = headshape;
cfg.coordsys = 'neuromag';
cfg.headshape.icp = 'no';

mri_realigned_3 = ft_volumerealign(cfg, mri_realigned_2);

% Save volume
save([output_path 'mri_realigned_3'], 'mri_realigned_3');

```

How does the head points align with the scalp now?

### Reslice MRI data
The final step of the co-registration is to reslice the MRI on to a 1x1x1 mm cubic grid aligned with the coordinate axes. 

```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.resolution = 1;

mri_resliced = ft_volumereslice(cfg, mri_realigned_3);
```

To avoid unit errors later in the pipeline, we also convert the units of the volume to centimeters. Then we save the volume.

```{r, engine='octave', eval=FALSE}
mri_resliced_cm = ft_convert_units(mri_resliced, 'cm');

save([output_path 'mri_resliced_cm'], 'mri_resliced_cm');
```

## Segment data into brain, skull, and scalp
We will not segment the volume into three compartments: One for the brain, one for the skull, and one for the scalp. This is done automatically with the function _ft_volumesegment_. All it needs is a configuration specifying which output segments we want.

```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.output = {'brain' 'skull' 'scalp'};

mri_segmented = ft_volumesegment(cfg, mri_resliced_cm);

save([output_path 'mri_segmented'], 'mri_segmented');
```

The following step is used to corcts the segmentation so we avoid overlap that will casue problems later on:

```{r, engine='octave', eval=FALSE}
binary_brain = mri_segmented.brain;
binary_skull = mri_segmented.skull | binary_brain;
binary_scalp = mri_segmented.scalp | binary_brain | binary_skull;
binary_scalp = mri_segmented.scalp + binary_skull;

% use boolean logic together with IMERODE
binary_skull = binary_skull & imerode(binary_scalp, strel_bol(2)); % fully contained inside eroded scalp
binary_brain = binary_brain & imerode(binary_skull, strel_bol(2)); % fully contained inside eroded skull
```

We will then copy the segmented volume into a new structure. We will also add the original anatomical image for plotting purposes.

```{r, engine='octave', eval=FALSE}
mri_segmented_2 = mri_segmented;  % Copy stucture
mri_segmented_2.anatomy = mri_resliced_cm.anatomy; % Coply anatomical image
  
% insert the updated binary volumes, taking out the center part for skull and scalp
mri_segmented_2.brain    = binary_brain;
mri_segmented_2.skull    = binary_skull & ~binary_brain;
mri_segmented_2.scalp    = binary_scalp & ~binary_brain & ~binary_skull;
```

Now it is a good idea to inspect that the segmentation went well. 

```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.funparameter = 'brain';
ft_sourceplot(cfg, mri_segmented_2);

cfg.funparameter = 'skull';
ft_sourceplot(cfg, mri_segmented_2);

cfg.funparameter = 'scalp';
ft_sourceplot(cfg, mri_segmented_2);
```

Brain segment:
![](./images/brain_seg.jpg)

Skull segment:

![](./images/skull_seg.jpg)
Scalp segment:

![](./images/skin_seg.jpg)

* Q: What is covered by the "brain" compartment?

## Construct meshes for each segment
From the three compartments we will now construct meshes that define the borders of the compartments; one for the outer scalp, one for the border between scalp and skull, and one for the inner skull, which we will call "brain".

To prepare the meshes, we use the FieldTrip function _ft_prepare_mesh_. We will also specify the number of vertices on each surface with the option _cfg.numverticies_.

```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.method = 'projectmesh';
cfg.tissue = 'brain';
cfg.numvertices = 3000;

mesh_brain = ft_prepare_mesh(cfg, mri_segmented_2);

cfg.tissue = 'skull';
cfg.numvertices = 2000;

mesh_skull = ft_prepare_mesh(cfg, mri_segmented_2);

cfg.tissue = 'scalp';
cfg.numvertices = 1000;

mesh_scalp = ft_prepare_mesh(cfg, mri_segmented_2);

% Collect meshes into a signle structure
meshes = [mesh_brain mesh_skull mesh_scalp];

% save
save([output_path 'meshes'], 'meshes');
```

Plot the meshes to see if they look fine:

![](./images/what_a_mesh_we_made.png)


### Construct separate head models for EEG and MEG
Now we are ready to construct our volume conduction model -- i.e. the head model. Head models define how electric potentials and magnetic fields spread throughout the conductor (the head). We need separate head models for source reconstruction of MEG and EEG data. 

For the EEG head model, we will use all three meshes to define a three-layer model consisting of three compartments. To do this, we specify the method as "singleshell" (_cfg.method = 'bemcp'_) for the function _ft_prepare_headmodel_ and give the structure containing all three meshes as input.

For MEG we will use a single volume model based on the "brain" mesh. To create the single volume model, we specify the method as "singleshell" (_cfg.method = 'singleshell'_) for the function _ft_prepare_headmodel_ and give only the brain mesh as input.

 
```{r, engine='octave', eval=FALSE}
cfg = [];
cfg.method = 'bemcp';
cfg.conductivity = [1 1/20 1] .* (1/3);

headmodel_eeg = ft_prepare_headmodel(cfg, meshes);

cfg = [];
cfg.method = 'singleshell';

headmodel_meg = ft_prepare_headmodel(cfg, mesh_brain);

save([output_path 'headmodel_eeg'], 'headmodel_eeg');
save([output_path 'headmodel_meg'], 'headmodel_meg');
```

Note that we define the conductivity of the compartments differently (the values in _cfg.conductivity:).

* Q: What could be a potential problem with how we define the conductivity of the different compartments, and how does the problem affect the MEG and EEG respectively?

Finally, we will plot the head models together with the head points and sensors. This is a sanity check that the process went well. Note the simplicity of the MEG model compared to the EEG model.

First plot the MEG:
```{r, engine='octave', eval=FALSE}
CODE;
```

![](./images/meg_check.jpg)

Then plot the EEG model. It can be difficult to see, but the figure should contain three meshes inside one another. If you rotate the figure, it might be easier to see.

```{r, engine='octave', eval=FALSE}
CODE;
```

![](./images/eeg_check.jpg)

