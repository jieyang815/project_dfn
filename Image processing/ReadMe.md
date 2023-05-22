# Image processing

## Data
Image data are saved under "Image" 
- Original image ("SingleImage....tif")
- Sample data used for reconstruction (5 regions of interests and corresponding ground truth)
- A simple macro for ImageJ used for manual segmentation ("MeasureAndGetCoor")
    

## Code

### MATLAB Toolboxes

All code is written in MATLAB. The following toolboxes are required:
- The Computer Vision Toolbox
- Image Processing Toolbox

### Main scripts

`main_optimize.m`: main algorithm performs fiber segmentation and parameter optimization based on bi-directional evaluation
`geoPreprocess.m`: preprocess the geometry of the reconstructed fiber network before discretization

### Self-defined functions

`FiberSeg.m`: segment fiber network from the images with preprocessing, fiber detection, and postprocessing.
`NetworkDetection.m`: detect nodes and trace fibers along centerlines using DFS.
`FiberConnect.m`: postprocess the fiber detection results to connect broken fibers.
`getMetric.m`: calculate the evaluation metrics based on Mayerich's method.
`myDetectHarrisFeatures`: self-defined function to detect nodes using Harris corner detection (makes some small changes based on the original in-built function).


## Procedure

### Reconstruct fiber network

1. Open `main_optimize.m`, and then run the algorithm.
2. Select the folder with images at the pop-up window (Images of 5 regions of interests are saved in "Sample data").
3. Then select the folder to save optimization results and the geometry of reconstructed fiber network.
4. The segmentation and evaluation will automatically begin afterwards. The whole process will take approximately 15 minutes for 5 images.
5. After the evaluation is completed, the optimal segmentation result for every image will be plotted, with the local evaluation metric mapped on the network.

(Note: the parameter optimization results are saved in "Optimization" for review. And the geometry of reconstructed fiber network is saved in "data_reNetwork")

### Preprocess the geometry of fiber network

1. Open `geoPreprocess.m`, and then run the algorithm.
2. Select the folder with the geometry of reconstructed network, and then select the folder to save the preprocessed geometry.
(Note: geometry of 3 networks is saved in the folder "Geometry_DFN", which can be used to test the algorithm)
3. The preprocessing will automatically begin afterwards. The whole process will take less than 1 second.

(Note: the preprocessed geometry is saved in "Geometry_DFN", which can be directly used to test the discretization algorithm in the "Numerical simulation" folder)
