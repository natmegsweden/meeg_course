# From raw data to evoked responses: pre-process MEG & EEG data

In this tutorial you will go though the various steps that takes you from raw data to evoked responses. These steps are:
1. Import the raw MEG and EEG data.
2. Inspect the data and metadata.
3. Select timelocked data
4. Clean the data by removing artifacts, filter data, and re-reference EEG.
5. Average trials to get evoked responses.

The raw data is stored in the three `fif` files that you have downloaded before beginning the tutorials:

    'tactile_stim_raw_tsss_mc.fif'
    'tactile_stim_raw_tsss_mc-1.fif'
    'tactile_stim_raw_tsss_mc-2.fif'
    
Note that though there are three data files, they are all part of one single recording session. But because the `fif` format (the file format of the Neuromag MEG/EEG system that we used to record data) only allows files to have a size of up to 2GB the recording has been split into separate files. We should therefore think of these as one single file when we continue to process the data despite there being three files.

This is an important thing to notice when recording data. It illustrates the need for consistent file naming. For the example data, you can see that all are called `tactile_stim`, indicating that it is the same task.

## Setup paths
The first step is to point to the path where we have the data and setup FieldTrip. Change these to appropriate paths for your operating system and setup

```matlab
clear all
close all
restoredefaultpath
addpath('C:/fieldtrip/')            % Change to match your FieldTrip path
ft_defaults

meg_path = 'C:/meeg_course/data';   % Change to match your data path
```

Then define subject and recording specific paths. For now we only have one subject and session. In principle, we could just define the path as one string variable when we only have one subject. But we introduce this already now as it is a good ´way to organize your data when you have multiple subjects and/or session. In that case, the cell array `subjects_and_dates` can be expanded to include more subjects, simply by adding the subject ids and session names.

```matlab
%% Define subject paths
% List of all subjects/session
subjects_and_dates = ...
                    {
                        'NatMEG_0177/170424/'  % add more as needed
                    };
           
% List of all filenames that we will import                
filenames = {...
        'tactile_stim_raw_tsss_mc.fif' 
        'tactile_stim_raw_tsss_mc-1.fif'
        'tactile_stim_raw_tsss_mc-2.fif'
            };

n_filenames = length(filenames);  

output_path = fullfile(meg_path, subjects_and_dates{1});
```

## First look at what is in the datafiles
The `fif` files contain everything that was recorded during the recording data, including MEG data, EEG data, triggers, and various metadata. Before you import everything, take a look at what is in the files. This is especially a good idea if you are dealing with large files to avoid that you confidentially read in more data that what your computer can handle.

Now it is time to use the first FieldTrip function: use `ft_read_header` to read metadata from the `fif` files. Note that this will not read the data yet. `ft_read_header` is a low-level FieldTrip function and does not need a `cfg` struct to work.

```Matlab
%% Read header
infile = fullfile(output_path, filenames{1});
hdr = ft_read_header(infile);
```
The output of `ft_read_header` is here the struct `hdr`. This contain information about what channels is in the data. Explore the struct to find out what is in the data file.

> **Question 1.1:** What type of data channels is in the data and what is the sample frequency? 

One thing to be aware of is that this only read the header information from one of the three `fif` files. The information about which channels are in the data is the same for all three files as it was recording in one session. But the duration and time-stamps of the recordings might be off because of the split files.

Now read the headers of all three data files to find out how long the recording was and how much data we have. Make a preallocated cell array and loop over data files to read the header information into the cell array:

```Matlab
% Read all headers
hdrs = cell(3,1);

for ii = 1:length(filenames)
    infile = fullfile(output_path, filenames{ii});
    hdrs{ii} = ft_read_header(infile);
end
```
Take a look at the information in the `hdrs` array (index cell arrays with curly brackets like this: `hdrs{1}`). Look at the channel information as you did before. The channel information should be the same for all three files.

> **Question 1.2:** how many samples are there in the total of the entire recording? How long is the time of the recording?
> 
>Hint: look at the fields `nSamples` and `Fs` to calculate the total duration of the recording session.

## Read trigger values
The data consists of tactile stimulation to all five fingers of the right hand. When each stimulation to a finger occurred is marked by a trigger in the data. We will use these triggers to select the parts of the data that we will analyse later on.

What the values of the triggers represents is something you always want to write down in a trigger manual or protocol so that you always know what the values represent. You can read the values from the data, but their meaning is something you should know. 

For this data, I know from my recording notes the trigger values represent the following:

    1  = Little finger tactile stimulation
    2  = Ring finger tactile stimulation
    4  = Middle finger tactile stimulation
    8  = Index finger tactile stimulation
    16 = Thumb tactile stimulation
    32 = New block begins
    64 = End of experiment
    
But knowing what the values represents is one thing. Another is to see how they actually look in the data. It is a good quality check to inspect how the trigger values appear in the data. For example, if we are to pilot a newly designed experiment, we want to make sure that the value and the order of the triggers appear correct.

To inspect triggervalues we use `ft_read_event`. Similar to `ft_read_header` this is also a low-level FieldTrip function and does not need a `cfg` struct to work. Now use `ft_read_event` to read the events (i.e. triggers) in the file you specified above:

```matlab
eve = ft_read_event(infile);
```
Look at the `eve` structure:

> **Question 1.3:** What are the values and the types of the events in `eve` and how many events are there in total?

Because there are several trigger channels in the data, FieldTrip read all channels as independent channels and combine everything into a single struct. There are advantages of this, but for now we are only interested in the composite trigger channel called `STI101`.

Now use the header info you read above to find the channel index of the channel called `STI101`:

```matlab
% Find only relevant channel
chanindx = find(~cellfun(@isempty, strfind(hdr.label, 'STI101')));

% Read events only from relevant channel
eve = ft_read_event(infile, 'chanindx', chanindx);
```
Again, be aware that this only read the events from one of the three `fif` files. The information about which channels are in the data is the same for all three files as it was recording in one session. As before, now read the events from all three data files to find all the triggers in the entire recording.

```matlab
% Read all events
eves = cell(3,1);
for ii = 1:length(filenames)
    infile = fullfile(output_path, filenames{ii});
    eves{ii} = ft_read_event(infile, 'chanindx', chanindx);
end

% Combine
eve = [eves{1}, eves{2}, eves{3}];
```
Now that we have all events combined in one struct, take a look at how many we have of each trial:

```matlab
% See unique events
unique([eve.value])

% Make a summary table of event count
[vals,~,idx] = unique([eve.value]);
n  = histc([eve.value], vals);
evetab = [vals; n];   % Make a quick table
```
> **Question 1.4:** What is the count of each trigger value?

In addition to knowing how many trials we have of each type, we also want to know how the trials are distributed over time. The time the trigger occurred is stroed in `eve.sample`.

If you, however, look at the minimum and maximum sample values in `eves{1}`, `eves{2}`, and `eves{3}`, e.g. like this:

```matlab
min([eves{1}.sample])   #Change "1" to "2" and "3"
max([eves{1}.sample])   #Change "1" to "2" and "3"
```
you will notice that the range of the samples is the same for all files. This would mean that the triggers for each part of the occurred within the same duration of time. This is obviously not possible. What has happened here is that the sample info is defined relative to the onset of the individual files rather than the onset of the recording session. Let us correct that before we plot the triggers across time:

```matlab
% Correct sample info for split files
sam1 = [eves{1}.sample]; 
sam2 = [eves{2}.sample]+hdrs{1}.nSamples;
sam3 = [eves{3}.sample]+hdrs{1}.nSamples+hdrs{2}.nSamples;
allsam = [sam1, sam2, sam3];
```
Then plot the trigger values across time:

```matlab
% Plot
figure
scatter(allsam, [eve.value])
```

Image


## Read raw data
Set preprocessing options

Set preprocessing options

```matlab
events = 8; % right index finger
% to do all five fingers set events as below:
% events = [1 2 4 8 16]; %% right little, ring, middle, index, thumb
```

* Neuromag data files (.fif) are split into several files due to a limitation of 2 GB for each files
* Event code is read in (8), but note that all can be read in (1, 2, 4, 8, 16)
* We anticipate the measured delay of 41 ms between the trigger and the actual stimulation
* MEG data and EEG data are handled separately since EEG data need to be referenced to a reference “channel” (here the average of all channels)
 *   Not much preprocessing is done to the data (we’ll do that later for the respective analyses (timelocked and time frequency representation (TFR)))
    The split files are appended to one another to ease data handling

```matlab
%% Read raw data
split_files_MEG = cell(1, n_filenames);
split_files_EEG = split_files_MEG;

for filename_index = 1:n_filenames

    filename = filenames{filename_index};
    full_path = fullfile(meg_path, subjects_and_dates{1}, filename);
 
    % define trials
    
    cfg = [];
    cfg.dataset = full_path;
    cfg.trialdef.prestim = 1.609; % seconds % adjusted with 41 ms
    cfg.trialdef.poststim = 1.691; % seconds
    cfg.trialdef.eventvalue = events;
    cfg.trialdef.eventtype = 'STI101';
    
    cfg = ft_definetrial(cfg);
    
    % preprocess
 
    cfg.demean = 'no';
    cfg.lpfilter = 'no';
    cfg.hpfilter = 'no';
    cfg.dftfilter = 'no';
    cfg.channel = 'MEG';

    split_files_MEG{filename_index} = ft_preprocessing(cfg);
    
    cfg.channel = 'EEG';
    cfg.reref = 'no';
    cfg.refchannel = 'all';
    
    split_files_EEG{filename_index} = ft_preprocessing(cfg);
         
end
```

```matlab
% append split files
cfg = [];
preprocessed_data_MEG = ft_appenddata(cfg, split_files_MEG{:});
preprocessed_data_EEG = ft_appenddata(cfg, split_files_EEG{:});
```

## Inspect raw data

Try to inspect the data and see if you can identify some bad electrodes (and MEG sensors)


cfg = [];
cfg.continuous = 'yes';
cfg.viewmode = 'vertical'; % also try 'butterfly'

ft_databrowser(cfg, preprocessed_data_MEG);

ft_databrowser(cfg, preprocessed_data_EEG);

