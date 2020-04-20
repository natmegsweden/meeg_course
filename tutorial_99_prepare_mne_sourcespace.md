# Prepare source space for MNE
  
## Prepare surface source space

Minimum-norm estimate (MNE) use a model of the cortical surface as source model. MNE can also be done in the toolbox MNE-C and MNE-Python. The pipeline for MNE was developed in collaboration with the people behind the MRI processing toolbox Freesurfer (http://freesurfer.net/). To create the source model to do MNEs we, therefore, need to export our MRI to a Freesurfer compatible format and then process the image with Freesurfer to get the surface model of the white-grey matter boundary. Then we use MNE-C (the toolbox) to sample the surface to a homogenous set of point.

To do this, you need to have Freesurfer and MNE-C installed on your PC. You can find out how to do this here (MNE-C): http://martinos.org/mne/dev/install_mne_c.html and here (Freesurfer) http://freesurfer.net/fswiki/DownloadAndInstall. However, this is not required for this workshop. The Freesurfer and MNE toolboxes are called directly from the terminal in Linux or Mac. Unfortunately, it cannot run on Windows.

Do not run the following on the workshop, but do read through to get an understanding of the processing steps.

## Export MRI to Freesurfer
First, we load the raw imported MRI from the first tutorial. As before we use ``ft_volumealign`` to convert the coordinate system of the MRI, but this time we convert it to the SPM/MNI coordinate system. Note that this coordinate system is defined differently than the *Neuromag* coordinate system we used previously. The origin of the MNI coordinate system is the anterior and posterior commissure.


```matlab
%% Prepare MRI for Freesurfer w/SPM
load('mri.mat') 

% Convert to spm coordinate system
cfg = [];
cfg.method    = 'interactive';
cfg.coordsys  = 'spm';
mri_spm    = ft_volumerealign(cfg, mri);

```

When the coordinate system is defined, reslice the image to obtain the resolution Freesurfer require.

```matlab
%% reslicing into FreeSurfer friendly coordinates
cfg = [];
cfg.resolution = 1;
cfg.dim = [256 256 256];

mri_resliced_spm = ft_volumereslice(cfg, mri_spm);

save('mri_resliced_spm.mat','mri_resliced_spm')
disp('Done')
```

We need the transformation matrix from raw voxels to spm coordinates in one of the following steps. We save this for later. However, this is a manually set, so if you save this yourself, you might get problems with source space alignment later on, unless you intend to run the whole procedure. You should find a file called "transform_vox2spm.mat" amongst the files you have downloaded.

```matlab
transform_vox2spm = mri_resliced_spm.transform;

save('transform_vox2spm.mat','transform_vox2spm');
```

Then use ``ft_volumewrite`` to save the MRI as a "mgz" file that Freesurfer will read. If you run this yourself, remember to change the path and filename. Here the filename is "sub02" with filetype "mgz".

```matlab
%% Write MRI
cfg = [];
cfg.filename = 'workshop_material/data/mri/freesurfer/Sub02';
cfg.filetype = 'mgz';
cfg.parameter = 'anatomy';

ft_volumewrite(cfg, mri_resliced_spm);
```

Freesurfer has a long a convoluted processing pipeline. We will help Freesurfer by segmenting the MRI using ft_volumesegment to create a brain mask and export it as we just have done. You can visualize the brain mask as before with ft_sourceplot.


```matlab
%% Save brainmask for Freesurfer
%Save mask
mri = ft_read_mri('workshop_material/data/mri/freesurfer/Sub02.mgz');
mri.coordsys = 'spm';

cfg = [];
cfg.output = 'brain';
seg = ft_volumesegment(cfg, mri);
mri.anatomy = mri.anatomy.*double(seg.brain);

cfg             = [];
cfg.filename    = 'workshop_material/data/mri/freesurfer/Sub02/sub02mask';
cfg.filetype    = 'mgz';
cfg.parameter   = 'anatomy';
ft_volumewrite(cfg, mri);
```

## Run Freesurfer pipeline

The next step is to process the data using Freesurfer. This step takes about 10 hours. Do not run this now! Instead just read through the next few steps. If you are familiar with Freesurfer, this should be familiar. If not, then just keep this as a reference.

Freesurfer is not a MATLAB toolbox, but a stand alone program. To process the MRI in Freesurfer, we have to open a Terminal. In the Terminal we set up the Freesurfer environment a:
  
  ```bash
# Set up freesurfer
export FREESURFER_HOME='/opt/freesurfer' #Where Freesurfer is installed
export SUBJECTS_DIR='/workshop_material/data/mri/freesurfer/Sub02' #Folder containing all subjects (output from Freesurfer)
export SUBJECT='Sub02'  # Name of the subject we will run now

source $FREESURFER_HOME/SetUpFreeSurfer.sh #setup Freesurfer
```

Freesurfer has a defined way of naming and organizing files. We, therefore, create a folder for our subject (here called "Sub02") that follows this convention. Then copy the MRI files "sub02.mgz" and "sub02mask.mgz" to the MRI folder in the subject folder and convert them using "mri_convert". Again, this is done from the Terminal.

```bash
mksubjdirs $SUBJECTS_DIR/$SUBJECT #Make a dir for our subject

# Copy files to folders
cp sub02.mgz $SUBJECTS_DIR/$SUBJECT/mri/sub02.mgz
cp sub02mask.mgz $SUBJECTS_DIR/$SUBJECT/mri/sub02mask.mgz

cd $SUBJECTS_DIR/$SUBJECT/mri

mri_convert -c -oc 0 0 0 sub02.mgz orig.mgz
mri_convert -c -oc 0 0 0 sub02mask.mgz brainmask.mgz

```

We then run a selection of the typical Freesurfer processing pipeline only omitting the first few steps, which we already did by importing the dicom files and creating a brain mask.

```bash
#Running selection of Freesurfer reconstruction (alternative to use recon-all -all)
recon-all -talairach      		-subjid $SUBJECT
recon-all -nuintensitycor 	-subjid $SUBJECT
recon-all -normalization  	-subjid $SUBJECT
recon-all -gcareg         		-subjid $SUBJECT
recon-all -canorm        		 -subjid $SUBJECT
recon-all -careg          		-subjid $SUBJECT
recon-all -careginv      		 -subjid $SUBJECT
recon-all -calabel       		 -subjid $SUBJECT
recon-all -normalization2 	-subjid $SUBJECT
recon-all -maskbfs       		-subjid $SUBJECT
recon-all -segmentation   	-subjid $SUBJECT
recon-all -fill           		-subjid $SUBJECT
recon-all -tessellate 		-subjid $SUBJECT
recon-all -smooth1    		-subjid $SUBJECT
recon-all -inflate1   		-subjid $SUBJECT
recon-all -qsphere    		-subjid $SUBJECT
recon-all -fix        		-subjid $SUBJECT
recon-all -white      		-subjid $SUBJECT
recon-all -finalsurfs 		-subjid $SUBJECT
recon-all -smooth2    		-subjid $SUBJECT
recon-all -inflate2   		-subjid $SUBJECT

cp $SUBJECTS_DIR/$SUBJECT/mri/sub02.mgz $SUBJECTS_DIR/$SUBJECT/mri/rawavg.mgz
recon-all -autorecon3 -subjid $SUBJECT

echo 'Done'
```

If you want to know more about what each step does see: http://surfer.nmr.mgh.harvard.edu/fswiki/ReconAllDevTable

## Setup source space with MNE-C
After running all of this, the Freesurfer subject folder will contain many different files. Most of this we do not need. What we need is the  file describing the whit-grey matter boundary surface. This file will be read by the MNE function ``mne_setup_source_space`` which will create a set of points that covers the cortical surface. Again, this is executed from the Terminal on a Linux or Mac PC. As with Freesurfer, we need to set up MNE-C before creating the source space:
  
```bash
#Now we will continue to make source space with MNE
export MNE_ROOT='/opt/MNE' #path to MNE
source $MNE_ROOT/bin/mne_setup #setup MNE

mne_setup_source_space --ico -6 #Create source space
echo 'Done'
```

This will create a file called ``Sub02-oct-6-src.fif`` containing points along the cortical surface. The ``--ico -6`` flag denotes the algorithm used for sampling the points. The one we used here gives us 4098 points from each hemisphere.

Now the Freesurfer/MNE-C part is done and we can return to MATLAB.
You will find the file ``Sub02-oct-6-src.fif`` in the tutorial datafiles.

## Create head model

Now we have our source model we need to create a volume model. We load the spm converted MRI from before. This is now in a different coordinate system than the sensors, so first we need to convert the coordinate system to the "neuromag" coordinate system using "ft_volumerealign". Note how the Neuromag coordinate system is defined when you do the realignment!

```matlab
load('mri_resliced_spm.mat')

cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'neuromag';
mri_resliced_neuromag = ft_volumerealign(cfg, mri_resliced_spm);
````

The next step is then the same as you did earlier: use ``ft_volumrealign`` to align MRI to head points measured in the MEG.

```matlab
input_path  = fullfile(meg_path, subjects_and_dates{1}, filenames{1});

headshape = ft_read_headshape(input_path);
cfg = [];
cfg.method = 'headshape';
cfg.headshape.headshape = headshape;
cfg.headshape.icp = 'yes';
cfg.coordsys = 'neuromag';

mri_realigned_neuromag2 = ft_volumerealign(cfg, mri_resliced_neuromag);

save([output_path 'mri_realigned_mne2'], 'mri_realigned_2');

cfg = [];
cfg.method = 'headshape';
cfg.headshape.headshape = headshape;
cfg.coordsys = 'neuromag';
cfg.headshape.icp = 'no';

mri_realigned_neuromag3 = ft_volumerealign(cfg, mri_realigned_neuromag2);

save('mri_realigned_mne3', 'mri_realigned_neuromag3');
````

When the MRI and head points have been aligned, you can create the head model as before, by segmenting the MRI into brain, skull and scalp compartments and create meshes of each.

```matlab
%% Segment MRI and make headmodel
cfg = [];
cfg.output = {'brain' 'skull' 'scalp'};

mri_segmented_mne = ft_volumesegment(cfg, mri_realigned_neuromag3);
mri_segmented_mne.anatomy = mri_realigned_neuromag3.anatomy;

save('mri_segmented_mne.mat', 'mri_segmented_mne');
disp('Done');

%% Apply correction to avoid overlap etc, in three-layer model
binary_brain = mri_segmented_mne.brain;
binary_skull = mri_segmented_mne.skull | binary_brain;
binary_scalp = mri_segmented_mne.scalp | binary_brain | binary_skull;
binary_scalp = mri_segmented_mne.scalp + binary_skull;

% use boolean logic together with IMERODE
binary_skull = binary_skull & imerode(binary_scalp, strel_bol(2)); % fully contained inside eroded scalp
binary_brain = binary_brain & imerode(binary_skull, strel_bol(2)); % fully contained inside eroded skull

mri_segmented_mne2 = mri_segmented_mne;
% insert the updated binary volumes, taking out the center part for skull and scalp
mri_segmented_mne2.brain    = binary_brain;
mri_segmented_mne2.skull    = binary_skull & ~binary_brain;
mri_segmented_mne2.scalp    = binary_scalp & ~binary_brain & ~binary_skull;

% Plot
cfg = [];
cfg.funparameter = 'brain';
ft_sourceplot(cfg, mri_segmented_mne2);
cfg.funparameter = 'skull';
ft_sourceplot(cfg, mri_segmented_mne2);
cfg.funparameter = 'scalp';
ft_sourceplot(cfg, mri_segmented_mne2);

%% construct mesh from brain
cfg = [];
cfg.method = 'projectmesh';
cfg.tissue = 'brain';
cfg.numvertices = 3000;
mesh_brain = ft_prepare_mesh(cfg, mri_segmented_mne);

cfg.tissue = 'skull';
cfg.numvertices = 2000;
mesh_skull = ft_prepare_mesh(cfg, mri_segmented_mne);

cfg.tissue = 'scalp';
cfg.numvertices = 1000;
mesh_scalp = ft_prepare_mesh(cfg, mri_segmented_mne);

meshes_mne = [mesh_brain mesh_skull mesh_scalp];

save('meshes_mne.mat', 'meshes_mne');
disp('Done');

% Plot meshes
figure
ft_plot_mesh(mesh_brain, 'edgecolor', 'none', 'facecolor', 'r')
ft_plot_mesh(mesh_skull, 'edgecolor', 'none', 'facecolor', 'g')
ft_plot_mesh(mesh_scalp, 'edgecolor', 'none', 'facecolor', 'b')
alpha 0.3
````

Finally, make a single shell model of the brain for MEG and a three-layer model for EEG.

```matlab
%% Make headmodel for MEG
% EEG
cfg = [];
cfg.method = 'bemcp';
cfg.conductivity = [1 1/20 1] .* (1/3);
headmodel_mne_eeg = ft_prepare_headmodel(cfg, meshes_mne);

% MEG
cfg = [];
cfg.method = 'singleshell';
headmodel_mne_meg = ft_prepare_headmodel(cfg, mesh_brain);

save headmodel_mne_eeg.mat headmodel_mne_eeg
save headmodel_mne_meg.mat headmodel_mne_meg

figure; ft_plot_vol(headmodel_mne_meg, 'facealpha', 0.5, 'edgecolor', 'none');
figure; ft_plot_vol(headmodel_mne_eeg, 'facealpha', 0.2, 'edgecolor', 'none','facecolor','b');
```