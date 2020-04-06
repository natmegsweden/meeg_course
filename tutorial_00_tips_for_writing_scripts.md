# Tips for writing effective analysis scripts

Writing data analysis scripts can easily become a mess! there are many steps to go from raw MEG/EEG data to the final results. It is important to keep track of what processing step that goes before and after another. Know what data that should be read in one step, saved to memory, and then read in the next step. If we mess this up, we might end up with invalid results. And it is easy to make errors, e.g., read wrong data files, having different versions of toolboxes, working in the wrong directory, etc., especially in MEG/EEG data processing where there are several manual steps and we often go back to re run analysis.

Below you find a quick list of recommendations to make it easier for you to write effective analysis pipelines. The recommendations are based on van Vliet (2019)[1] and the MEG-BIDS guidelines[2]. I recommend that you take a look at these when you have to write your own analysis scripts.

## Use section breaks when testing code
When writing code you often want to run only a small snip of code, e.g. when you test your code while scripting.

If you only want to run parts of your script, you can mark code and select `run selection` or press `F9`. However, constantly marking code manually becomes annoying really fast. Instead, use section breaks. You start a section with `%%` (two percentage signs). The line is commented and you will notice that the section is highlighted in a yellowish color. If you now press `ctrl+enter` you will run the code in the highlighted section.

```Matlab
%% Make a section
x = 1:10;
y = rand(size(x));

%% Make a new section
plot(x, y, 'or')
```

## Define the paths and toolboxes in the beginning of the script
For this tutorial you will use the toolbox FieldTrip to analyses MEG/EEG data. FieldTrip is written in MATLAB but is not a part of MATLAB. We therefore need to make sure that MATLAB has FieldTrip in its PATH definition to use the functions. The same apply if you use other toolboxes. This is simple: simply use the MATLAB function  `addpath( ... )` to add the path where you downloaded FieldTrip. If you have several versions of FieldTrip or have used other toolboxes before you run this script, it is also a good idea to restore the PATH with the function `resotredefaultpath`

The start of my script my look like this:

```Matlab
resotredefaultpath
addpath('/home/mikkel/fieldtrip')
ft_defaults
```

If you have several variables with the same names, it might also be good to add 
```matlab
close all
clear all
````
before the above code, to clear all figures and variables from your workspace.

After this, we are ready to begin our script and we will use the version of FieldTrip that we know we have at the given location.

## Run all analysis with one version of the software
Toolboxes for data analysis gets updated regularly. FieldTrip, for example, get a new version daily. This is to keep the functionality up to date or to fix bugs in the code that users might have occurred. However, do not update FieldTrip daily! When you begin a project, make sure that all data is processed with the same version of the software that you use. If you need to update, which you sometime need to do, make sure that your code is backwards compatible with the updated toolbox. If you need to update, better run everything again.

This is why it is useful to have several versions to make sure that you use the same versions for each project.

## One script does one data processing step
It might seem like a good idea to have one big script that you only have to run once to go from raw data to finished result. It is not! It only makes it difficult to find bugs and errors. Instead try to follow the principle:

> **one script = one analysis step**

For example, one script that import raw data from disc, does the pre-processing, and then save the processed data to script. You can then easily name your scripts in the order they should run and call them from master script.

## Save intermediate results
Save data after each processing step. Especially before and after steps that require manual intervention. This makes it easier to go back and redo individual steps without the need to re-run everything again.

## Visualize results of intermediate processing steps
Thought the tutorial you will be asked a lot to plot data and inspect data structures. This is not just to keep you occupied! It is good practice to visuale the outcome of each processing step (when you also would save the data). When you visualize the output of each processing step, it is easy to spot errors as they occur rahter than only noticing them in the end results—if you even notice them at all by then.

## Use consistent filenames
Do not rename files each time you run the analysis. Use a consistent way to easy read what subject, session, and processing step the data belongs to. For example, output files from an analysis might look like this:

    sub01-raw-tsss.fif
    sub01-raw-downsampled.mat
    sub01-epochs.mat
    sub01-tfr.mat
    ...
    
Note that each file has the id of the subject (`sub01`) in all filenames and a string indicating what analysis step it belongs to. 

## Store data separate by subject and session
When you have data from multiple subjects resist the temptation to throw all data into one folder. Instead, create a project folder where you have one folder per subjects. And if you have more than one session per subject, you should then have separate sub-folders in the subject folder:

    /home/mikkel/my_project/data/ ...
        ./sub01 ...
            ./session1 ...
                ./sub01-ses1-raw-tsss.fif
                ./sub01-ses1-raw-downsampled.mat
                ./sub01-ses1-epochs.mat
                ./sub01-ses1-tfr.mat
                ... etc.
            ./session2 ...
                ./sub01-ses2-raw-tsss.fif
                ./sub01-ses2-raw-downsampled.mat
                ./sub01-ses2-epochs.mat
                ./sub01-ses2-tfr.mat
                ... etc.
        ./sub02 ...
        ... etc.

When you keep this structure, it is  easy to set up subject specific paths to loop though when processing multiple subjects.

In the tutorial data, you will find one subject called *NatMEG_0177* with one session called *170424* (the recording date; not the best session name). 

We can then setup subject and recording specific paths as below. The cell array `subjects_and_dates` can be expanded with more subjects and sessions when needed. You will see these lines of code several times throughout the tutorials.

```Matlab
%% Define subjects and sessions

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

In your scripts, you can then easily loop though several subjects and run the same processing step on all subjects. You also make sure that you always read and save data to the correct folder by generating the paths within the loops, rather than specifying it manually in each script; e.g. like this:

````Matlab
output_path = fullfile(meg_path, subjects_and_dates{1});
````

## Specify paths and subject names once
When you run multiple scritpts on several subjects and have to go back to redo steps of the analysis you can easily loose track of what files belong to which and when. You might end up with several versions of the same file processed with slightly different options. You want to make sure that when you go on to the next step, you read the correct files. If you change the filenames along the way, it is easy to loose track of which files that are new and old and you are prone to make ciritical errors

The best way to avoid such errors is to specifiy subject id and filenames as few places as possible; ideally only once! This can be done by making a meta-script where you specify filenames and subject id (such as the code snippet above) that you then run in the beginning of each script.


***

[1]: 

> van Vliet, M. (2019). [Guidelines for data analysis scripts](https://arxiv.org/pdf/1904.06163.pdf). *ArXiv:1904.06163*. 

[2]: [The Brain Imaging Data Structure (BIDS)](https://bids.neuroimaging.io/) is an initiative to standardise how neuroimaging data is stored  to allow easy sharing of data across reserach groups and sites. Originally BIDS was for sharing MRI data but has since been expanded for MEG (and by extension EEG). 
    
> Niso, G., Gorgolewski, K. J., Bock, E., Brooks, T. L., Flandin, G., Gramfort, A., Henson, R. N., Jas, M., Litvak, V., T. Moreau, J., Oostenveld, R., Schoffelen, J.-M., Tadel, F., Wexler, J., & Baillet, S. (2018). [MEG-BIDS, the brain imaging data structure extended to magnetoencephalography](https://doi.org/10.1038/sdata.2018.110). *Scientific Data, 5(1)*. 
