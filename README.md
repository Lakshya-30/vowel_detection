# vowel_detection
AIM: Vowel Recognition Assignment
EXECUTION: Build and run the code on Visual Studio 2010. Use F5 key to run the code.

INPUT: Training and testing speech samples can be found in the vowels directory and can be accessed by using 
	fopen_s(&f1, "vowels/<filename>.txt", "r") during file reading.
OUTPUT: The AvgCi's are stored in AvgCi directory and the predictions & accuracy are printed on the command prompt.

CONSTANTS:
1. framesize: The size of each frame (320 samples)
2. P: The number of past samples to consider (12)
3. F: Number of stable frames to consider (5)
4. tokhuraWeights: Array to store the 12 Tokhura weight values.
5. vowels: Array to store the 5 vowels.

VARIABLES:
1. samples: Array to store input samples.
2. energy: Array to store computed energy for each frame.
3. steadyFrames: Array to store the steady frames.
4. R: Array to store the Ri values.
5. A: Array to store the Ai values.
6. limit: The max amplitude value for normalization.
7. dcShift: DC offset.
8. nFactor: Normalization scale.
9. sampleSize: Total number of samples.
10. steadyStart: Starting sample of the steady frames.
11. steadyEnd: Ending sample of the steady frames.
12. C: Array to store the Ci values of 1 frame
13. avgCi: Array to store the frame wise average ci values of all vowels.
14. allCi: Array to store 100 rows of ci values.
15. refCi: Array to store the reference Ci values of the vowels.
16. tokhuraDist: Array to store the tokhura distance values.
17. accuracy: Array to store the accuracy values for each vowel.


FUNCTIONS:
1. getDCShift():
	Function to compute the DC shift (mean value) of the audio samples from the input file.

2. getNFactor():
	Function to compute the normalization scaling factor.

3. findSteadyFrames():
	Function to identify stable frames from the audio samples based on energy.

4. hammingWindow(): 
	Function to apply a hamming window to the frame-wise samples to reduce spectral leakage in the frequency domain.

5. computeRis():
	Function to compute the autocorrelation values R[f][m] for each frame.

6. computeAis():
	Function to calculate the LPC coefficients from the autocorrelation values using Levinson-Durbin algorithm.

7. raisedSinWindow():
	Function to apply the raised Sin Window to the Ci's of each frame

8. storeCito3D():
	Function to store the values of Ci for 1 frame in the allCi array

9. computeAvgCi():
	Function to store the avg ci values to file for reference.

10. training():
	Function to read the first 20 instances of each vowel and perform training.

11. tokhuraDistance():
	Function to calculate the Tokhura distance using dump Ci values of training set.

12. calculateTokhura():
	Function to make prediction based on the Tokhura distance

13. testing():
	Testing the last 10 utterances of each vowel.


PROCEDURE:
Generating the reference file for a vowel:
1. Take the first vowel recording, and perform a DC shift and normalization.
2. Select 5 frames from the steady part.
3. Compute the Ri's, Ai's, and Ci's for these frames, and apply a raised Sine window to the Ci's.
4. Repeat these steps for 19 more recordings of the same vowel.
5. This will give you 100 rows of Ci values for the vowel:
(20 recordings   5 frames per recording   12 Ci values per frame) = 100 rows.
6. Average the Ci values for corresponding frames across all 20 recordings (e.g., average all first frames, all second frames, etc.).
7. This results in 5 rows of Ci values (5 rows   12 columns).
8. Save these values in a text file.
Repeat the process for the remaining vowels, resulting in 5 text files one for each vowel.

For testing:
1. Take 10 test files per vowel, and process them in a loop to check how many are recognized correctly.
2. For each test file take 5 frames from the stable part.
3. Calculate Tokhura's distance for each test file using the reference files.
4. Since each test file has 5 frames and each reference file has 5 rows, calculate Tokhura s distance for corresponding frames. Average the 5 distances to get the final distance.
5. Repeat this for each reference file.
6. The vowel with the smallest distance will be recognized as the correct vowel.
