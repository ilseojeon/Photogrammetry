## Go to master branch!

# Photogrammetry

---

### What are these codes about?

These are from a class "Advanced Photogrammetry" lectured by Prof. Impyeong Lee from LSM(Lab for sensor and modeling), dept. of Geoinformatics, University of Seoul.

I'm a graduate student of Geoinformatics and took this class in the spring semester, 2020.

I want to share my assignments for who's willing to start to learn basic codes of SfM(Structure-from-Motion), widely known as a photogrammetric process. The codes are written in MATLAB.

Instructions are below for each folder.

---

### Instructions for each folder

Lectures are organized in 7 labs, aiming to understand bundle adjustment which plays an important role in SfM frameworks(also used in SLAM nowadays). Each lab has the following main codes, showing the goal of it.

**Lab1**

Lab1_2a.m: Import and display drone images with EXIF tags (absolute positions and euler angles included).

Lab1_2b.m: Visualize camera coordinate systems of each images

**Lab2**

Lab2_model1.m: Visualize the ground coverage of images with photogrammetric model(the collinearity equation)

Lab2_model2.m: Visualize the ground coverage of images with computer vision model(camera projection matrix)

**Lab3**

Lab3_Q1.m: Develop a function of least squares estimation based on Gauss Markov Model(GMM) and drawing a straight fitting line with the function

Lab3_Q2.m: Develop the observation equations for the horizontal network adjustment to determine the locations of multiple cars

**Lab4**

Lab4_main.m: Determine a position of a single photo with ground control positions (called "Single Photo Resection")

**Lab5**

Lab5_main.m: Determine relative orientations with drone images from manually detected features

**Lab6**

lab6_fun.m: Estimate the fundamental matrix with matched images and assess the result

lab6_ho.m: Estimate the homography matrix with matched images and assess the result

**Lab7**

AT_Estimate_Ke.m: Implement bundle adjustment with only images

AT_Estimate_KeKg.m: Implement bundle adjustment with both images and ground control points

AT_Estimate_Kg.m:Implement bundle adjustment with only ground control points

Although each folder has the detailed information by .pdf files, those are written in Korean. I tried to make all codes clear with comments in codes instead. 

---

### Referenced by

"Advanced Photogrammetry" lectured by Prof. Impyeong Lee from LSM(Lab for sensor and modeling), dept. of Geoinformatics, University of Seoul.
