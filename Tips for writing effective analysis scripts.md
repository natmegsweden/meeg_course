# Tips for writing effective analysis scripts

Writing data analysis scripts can easily become a mess! there are many steps to go from raw MEG/EEG data to the final results. It is important to keep track of what processing step that goes before and after another. Know what data that should be read in one step, saved to memory, and then read in the next step. If we mess this up, we might end up with invalid results. And it is easy to make errors, e.g., read wrong data files, having different versions of toolboxes, working in the wrong directory, etc., especially in MEG/EEG data processing where there are several manual steps and we often go back to re run analysis.

Below you find a quick list of recommendations to make it easier for you to write effective analysis pipelines.

## Define the paths and toolboxes in the beginning of the script
For this tutorial you will use the toolbox FieldTrip to analyses MEG/EEG data. FieldTrip is written in MATLAB but is not a part of MATLAB. We therefore need to make sure that MATLAB has FieldTrip in its PATH definition to use the functions. The same apply if you use other toolboxes. This is simple: simply use the MATLAB function  `addpath( ... )` to add the path where you downloaded FieldTrip. If you have several versions of FieldTrip or have used other toolboxes before you run this script, it is also a good idea to restore the PATH with the function `resotredefaultpath`

The start of my script my look like this:

```octave
resotredefaultpath
addpath('/home/mikkel/fieldtrip/fieldtrip')
ft_defaults
```
After this, we are ready to begin our script and we will use the version of FieldTrip that we know we have at the given location.

## Run all analysis with one version of software
Toolboxes for data analysis gets updated regularly. FieldTrip, for example, get a new version daily. This is to keep the functionality up to date or to fix bugs in the code that users might have occurred. However, do not update FieldTrip daily! When you begin a project, make sure that all data is processed with the same version of the software that you use. If you need to update, which you sometime need to do, make sure that your code is backwards compatible with the updated toolbox. If you need to update, better run everything again.

This is why it is useful to have several versions to make sure that you use the same versions for each project.

## One script does one data processing step
It might seem like a good idea to have one big script that you only have to run once to go from raw data to finished result. It is not! It only makes it difficult to find bugs and errors. Instead try to follow the principle that one script = one analysis step. For example, one script that import raw data from disc, does the pre-processing, and then save the processed data to script. You can then easily name your scripts in the order they have to be run and call them from  master script.

## Save intermediate results
Save data after each processing step. Especially before and after steps that require manual intervention. This makes it easier to go back and redo individual steps without the need to re-run everything again.

## Use consistent filenames
Do not rename files each time you run the analysis. use a consistent way to easy read what subject, session, and processing step the data belongs to. For example, output files from an analysis might look like this:
```
sub01-raw-tsss.fif
sub01-raw-downsampled.mat
sub01-epochs.mat
sub01-tfr.mat
...
```
Note that each file has the id of the subject (`sub01`) in all filenames and a string indicating what analysis step it belongs to. 

## Store data separate by subject and session
When you have data from multiple subjects resist the temptation to throw all data into one folder. Instead, create a project folder where you have one folder per subjects. And if you have more than one session per subject, you should then have separate sub-folders in the subject folder:

````
/home/mikkel/my_project/data/ ...
    ./sub01 ...
        ./session1 ...
            ./sub01-ses1-raw-tsss.fif
            ./sub01-ses1-raw-downsampled.mat
            ./sub01-ses1-epochs.mat
            ./sub01-ses1-tfr.mat
            ...
        ./session2 ...
            ./sub01-ses2-raw-tsss.fif
            ./sub01-ses2-raw-downsampled.mat
            ./sub01-ses2-epochs.mat
            ./sub01-ses2-tfr.mat
            ...
    ./sub02 ...
    ...
````

When you keep this stucture, it is  easy to set up subject specific paths to loop though when processing multiple subjects. In the example below, we have a subject called *NatMEG_0177* with one session called *170424* (the recording date; not the best session name). 

We can then setup subject and recording specific paths as below. The cell array `subjects_and_dates` can be expanded with more subjects and sessions when needed. 

```octave
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
```

````matlab
output_path = fullfile(meg_path, subjects_and_dates{1});
````


## Quick tricks for using MATLAB

If you only want to run parts of your script, you can mark code and select "run selection" or press F9. However, constantly marking code manually becomes annoying really fast. Instead, use section breaks. You start a section with `%%`. The line is commented and you will notice that the section is highlighted in a yellowish color. If you now press `ctrl+enter` you will run the code in the highlighted section.
```matlab
%% Make a section
x = 1:10;
y = rand(size(x));

%% Make a new section
plot(x, y, 'or')
```