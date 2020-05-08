# HTDCountingProject
An attempt at using the simple python image processing pyplot module to count cells from raw microscope .tif data. It is similar to watershed, but uses my own unique algorithm in an attempt to get more accuracy.

DISCLAIMER: YES, I AM REINVENTING THE WHEEL HERE. THIS HAS BEEN DONE IN MUCH MORE EFFECTIVE WAYS BY PROGRAMMERS MUCH MORE EXPERIENCED THAN I AM, I AM ONLY A STUDENT LEARNING THE PYTHON LANGUAGE.

Project goal(s) and general overview: 

This program has the goal to perform a semiautomatic count of the number of cells, as well as how many cells (defined by blue nuclei) are marked as green and/or red (or none) using a combination of user input, with their human intuition being used to determine filter levels and maximum distance to make assignments, and Python to automate counting. -user manually adjusts the filter level of image to a level they feel most comfortable with.

Program will: 

-accept 3 black and white "channels" from the same microscopy image of cells" collected as raw tiff images, one of the blue channel, one red channel, one green channel. 
-count number of cell nuclei in an image with b/w tiff image of blue DAPI channel, parameters filtering for brightness beaks for counting will be defined by the user. We could employ special tactics to distinguish area sizes to get a more accurate cell count, but there are several mega cancer cells in my photo that should only be counted as once cell 
-mark down the coordinates of the center of the nuclei enclosed by smallest rectangles 
-count number of "blobs" of pigment using same method for cell nuclei, but instead the center of all blobs is assigned to the nearest nucleus. If the distance is larger than the parameter defined by the user, the blob will skip being assigned. 
-if color of blob is already assigned to the nucleus, it is only ever counted once. 
-cell nuclei coordinates are categorized based in whether they match at least one red blob, green blob, or both. 
-number of "none" marked nuclei are determined by subtracting total counted pigmented cells minus number of total nuclei. 
-totals are counted up and displayed to the user.

How to use the Program: 

1. Before running the python script, run each of the channel images provided for filtering, and write down the value that feels best to you. Make sure to write down the file names for each channel. Make sure the channel images are in the same folder as the python program. 

2. Look at the full image with channels combined, pick out the farthest cell marker from it's nuclei that you would personally count as being associated(not a fluke in staining) write down the approximate distance in pixels. 

3. Run the program 'ProjectHTD.py' from the command line by first redirecting to the command directory using command cd. 

4. Direct the command line to the location of the python program by pasting in the file directory of the program. 

5. Run the program by using the command "python ProjectHTD.py" 

6. Input the name of the files, the filtering values, and maximum distance when prompted. 

7. The output will be delivered to the folder as a txt file, or you can write down the values delivered from the command prompt.
