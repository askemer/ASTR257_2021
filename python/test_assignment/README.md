
# Overview

Both of your Fall quarter grad classes (ASTR 257 and ASTR 204) require Python programming.  If you feel that you would benefit from going through some brief Python tutorials, we have compiled some resources below.  To be clear, we expect many of you are proficient in Python (probably more proficient than your instructors), but if you’ve mostly programmed in other languages, it would be worthwhile to practice Python before classes start.

To define “proficient” we have a brief problem for you to try.  This is optional and 
you don’t have to turn it in.  If this is straight-forward for you, then you should have 
no problem with the programming for ASTR 257.

 
 ## Test assignment

1. Open the file test.fits in Python.  
2. Crop a 100x100 pixel region around the brightest star in the image.  
3. Save the cropped region as test2.fits.
4. Use the code *photutils* to measure the centroid of this star in the uncropped region

If you can’t do this problem, or would just like some more practice, 
we’ve compiled some tutorials below.
# Basic Python for Astronomy

There is a lab class at Amherst College that begins with an excellent set of Python tutorials:

[Amherst tutorials](https://github.com/spacegal-spiff/AST337-Fall2018)


In the folder AST337-Fall2018/Homeworks/JupyterExercises/ there are two exercises that take you through installing Python, programming in the Jupyter environment, basic mathematical programming, arrays and plotting.  This is an excellent place to start. 


Once you’ve done that, Lab 1 (/AST337-Fall2018/Labs/Lab01/) goes through an exercise to read in tabular data and plot it.

 
# I/O with .fits files


The same lab class has a lab that involves reading in a .fits file (standard astronomical image file) and displaying it (/AST337-Fall2018/Labs/Lab04/).

 

More generally, here is astropy’s description of its fits utility:

 
[astropy FITS/IO](https://docs.astropy.org/en/stable/io/fits/)