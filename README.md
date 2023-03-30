# AuSoMS - Automatic Segmentation of the Mouse Skull in MR images  
<img src="https://github.com/Weitspringer/ausoms/blob/main/thesis/Figures/04_exemp_seg.png" alt="Segmented mouse skull in MRI sequence" width="50%"/>

## Prerequesities
You will need a folder `<dicomFolder>` containing a DICOM sequence. AuSoMS was developed for slices in the transversal plane.
## Installation
Clone the project and open the project directory in MATLAB.
## Launch
Add the whole project folder `./segmentation` in the MATLAB path. The main code resides in `segmentation/functions`. \
To call the default segmentation, execute
```
[skull_segmentation, brain_segmentation, ~] = segmentVolume(<path-to-dicomFolder>, struct());
```
AuSoMS will return a skull segmentation as well as a brain segmentation. They are represented by a 3D matrix each. 

## Options
In order to run AuSoMS with custom options, create an empty struct 
```
options = struct();
```
and define attributes for it:
- `options.resizeFactor` Resizing factor of volume. Default: 1
- `options.SEBg` Structuring element used for background removal (**Caution: Structuring element will be used as is, ignoring the resize factor**)
- `options.SEBgClose` Structuring element used for closing of foreground (**Caution: Structuring element will be used as is, ignoring the resize factor**)
- `options.sizeSEBrain` Structuring element size for 3D PCNN brain segmentation in pixels (**Caution: Size will be used as is, ignoring the resize factor**)
- `options.cutoff` Number of slices ignored at the start as well as the end of the sequence. Default: 0
- `options.nSegmentsForeground` Number of segments for the foreground. Default: 4
- `options.brainSegmentation` Logical volume of existing brain segmentation (same size as input volume!). Alternative to using built in 3D-PCNN segmentation of the brain.
- `options.maxSkullThickness` Define thickness of skull (in mm). Default: 0.7
- `options.evalFile` _.mat_ file which contains ground truth of skull `skulltruth`, scale factor `factor` in relation to scale of input slices and, optionally, ground truth of the brain `braintruth`
- `options.fScoreBeta` Indicates how 'many times' recall is more important than precision. Default: 2
- `options.pcnnRepetitions` Repetitions for averaging of brain segmentation. Default: 1
- `options.outFile` Specify file for saving the results. Default: No file specified

You don't have to define every attribute. Defining an attribute solely overwrites parameters in the main function later on.
Then call
```
[skull_segmentation, brain_segmentation, ~] = segmentVolume(<path-to-dicomFolder>, options);
```
For example, if you want to work with the assumption that the maximum skull thickness is 0.6 mm instead of 0.7 mm (default), the code should look like this
```
options = struct();
options.maxSkullThickness = 0.6;
[skull_segmentation, brain_segmentation, ~] = segmentVolume(<path-to-dicomFolder>, options);
```
You can also save the results in a _.mat_ file by defining
```
options.outFile = <path-to-saveFile>;
```
That's it! For a more detailed explanation on these parameters, read the thesis in the `thesis` project directory.

## Evaluation
In case you want to evaluate the results, make sure you have a _.mat_ file `<evalFile>` containing the variables `skulltruth` as well as the rescale factor `factor` (use value 1 if image was not rescaled for manual segmentation) and - optionally - `braintruth`. Now execute
```
options = struct(); options.evalFile = <path-to-evalFile>;
[skull, brain, eval] = segmentationVolume(<path-to-dicomFolder>, options);
```
The variable `eval` contains the Dice Score and the Matthews Correlation Coefficient for the skull segmentation (whole as well as upper half) and the brain segmentation if a ground truth was given. The scores will also be saved in `<outFile>` if you defined `options.outFile` before executing the main function.
