# Discussion of making JRCLUST work with NDI

There are 2 ways to create the user experience. First, one might have the user launch JRCLUST and just
specify NDI data as the target of analysis. Second, one might have NDI launch JRCLUST.

Let's do the first one at this time.

Actions:

## Need to make an NDI recording object in JRCLUST (subclass of RawRecording)

Describe here

## Need to allow parameters to specify an NDI record:

New interpretations and variables needed for *.prm file:

```
		rawRecordings = {epoch IDs to include}
		recordingFormat = 'NDI';

		New variable: NDI (a structure)
		NDI.experimentpath = path to experiment;
		NDI.element_name = name of the element
		NDI.reference = reference number of the element

```

The following parameters need to be filled out:

```
	nChans = 32; % Number of channels stored in recording file (Distinct from the number of AP sites)
	sampleRate = 20000; % (formerly sRateHz) Sampling rate (in Hz) of raw recording

	probePad = [12, 12]; % (formerly vrSiteHW) Recording contact pad size (in μm) (Height x width)
	shankMap = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1]; % (formerly viShank_site) Shank ID of each site
	siteLoc = [0, 0; 0, 5; 0, 10; 0, 15; 0, 20; 0, 25; 0, 30; 0, 35; 0, 40; 0, 45; 0, 50; 0, 55; 0, 60; 0, 65; 0, 70; 0, 75; 0, 80; 0, 85; 0, 90; 0, 95; 0, 100; 0, 105; 0, 110; 0, 115; 0, 120; 0, 125; 0, 130; 0, 135; 0, 140; 0, 145; 0, 150; 0, 155]; % (formerly mrSiteXY) Site locations (in μm) (x values in the first column, y values in the second column)
	siteMap = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]; % (formerly viSite2Chan) Map of channel index to site ID (The mapping siteMap(i) = j corresponds to the statement 'site i is stored as channel j in the recording')

```

The following parameters can be "bootstrapped" (JRCLUST term for assigning default parameters) as follows:

		PREPROCESSING parameters can go through as defaults

		SPIKE DETECTION PARAMETERS can go through as defaults
			(User might want to edit nSiteDir)

		The following can all take their default values:
			FEATURE EXTRACTION PARAMETERS
			DISPLAY PARAMETERS
			AUX CHANNEL PARAMETERS
			LFP PARAMETERS
			TRACES PARAMETERS
			PREVIEW PARAMETERS
			VALIDATION PARAMETERS
			TRIAL PARAMETERS
			CLUSTERING PARAMETERS
			CURATION PARAMETERS

		

