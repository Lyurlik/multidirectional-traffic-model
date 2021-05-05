# Multidirectional Traffic Model

This is the project developed in scope of PhD thesis titled "Traffic control in large-scale urban networks" by Liudmila Tumash in DANCE team of GIPSA-Lab, CNRS, Grenoble, France. 

This code is used to produce two different vehicle density distributions. One density is the result of numerical simulation of NSWE model that is a macroscopic two-dimensional first-order PDE system consisting of four equations that predicts vehicle density evolution in four cardinal directions (North, South, West and East), and the other density is the one reconstructed from data obtained from real sensors and floating-car data on 8.01.2021 in Grenoble, France.

The project was built with Qt 5.12 MSVC 2017 64bit. In addition to Qt it requires Eigen library, path to which should be specified in Traffic_2D.pro file.

# Data files

In order to run the code, you need to have the following files:

1. Network's topology:

  * _"../ModelValidation/IntersectionTable.csv"_ -- contains information about intersections: x and y coordinates of every intersection (columns 1 and 2), its ID (column 3) and whether it is a node on border (column 4), which means that this intersection is located at domain's boundary through which vehicles may enter (inflows), or exit (outflows);

  * _"../ModelValidation/RoadTable.csv"_ -- contains information about roads: ID1 and ID2 (columns 3 and 4) are the id's of corresponding intersections that the road is connecting, ID_road (column 5) is the road's ID, max_vel (column 6) is its free-flow limit estimated from real measurements, then we have number of lanes (column 7) and road's length (column 8);

  * _"../ModelValidation/TurnTable.csv"_ -- contains turning ratios between any pair of roads: ID1 of incoming road (column 1), ID2 of outgoing road (column 2) and the turning ratio between these roads (column 5). 

2. Data from real sensors (model validation part):

  * _"../ModelValidation/Timestamp.csv"_ -- contains time in seconds at which the data are given (unix timestamp), the time step equals to one minute;

  * _"../ModelValidation/Density.csv"_ -- contains estimated density from real sensors: first number is road_id followed by its density (that is assumed to be constant within one road) at all time instants, then the next road_id with its density data for each time instant and so on;

  * _"../ModelValidation/AllInflows.csv"_ -- contains inflow values (in veh/hour) for every road for every time step (one minute). If road is outgoing from intersection that is not on border, then the inflow value is zero;

  * _"../ModelValidation/AllOutflows.csv"_ -- contains outflow values (in veh/hour) for every road for every time step (one minute). If road is incoming into intersection that is not on border, then the outflow value is zero.

 The data provided in this repository come from real measurements taken on 08.01.2021 in the city of Grenoble, France. The path to the files is specified in mainwindow.cpp.

# Code structure

The main file of the project is mainwindow.cpp: in its constructor we specify the file names to be loaded, start simulation starting time (line 26) and simulation step size (line 28). We can also change there the weighting parameter eta used to approximate parameters for every cell (line 4), and parameter d0 (line 5) is used for Gaussian Kernel estimation. 

Other important classes are:
* __UrbanNetwork__, which contains all the network geometry information (this is where all the network files are read). This network is used for both densities. In its function loadRoads one needs to specify the minimum distance between the heads of two consequative vehicles.
* __NSWEmodel__, which contains translation procedure of all network and intersection parameters into NSWE-formulation (function processIntersections). After all parameters are defined in NSWE, it calls constructInterpolation function that approximates these parameters defined for every intersection to be defined on every cell of a network. Then update is performed, where the Godunov numerical scheme is applied for the state update using NSWE model.
There is also a function getSSIMDiff_mean_weighted used to compute the weighted SSIM index between two densities (one from NSWE and the other one from real data).
* __GrenobleData__, where all the data estimated from the real-life experiments are loaded. In function reconstructDensity the density initially given for each road is defined for every cell. Thereby, every road is divided in 10 parts and density values are presented as points on the border between these parts. Then Gaussian Kernel estimation is used to determine density for every cell in the domain.
* __TrafficSystem__ implements concurrent thread for parallel NSWE simulation relative to the main visualization thread
