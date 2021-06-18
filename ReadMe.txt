This project is used to produce two different vehicle density distributions. One density is the result of numerical simulation of NSWE model that is a two-dimensional first-order PDE system consisting of four equations, and the other density is the one reconstructed from data obtained from sensors and floating-car data in Grenoble.

In order to run the code, you need to have the following files:

I. Read network's topology:

1) "../ModelValidation/IntersectionTable.csv" -- contains information about intersections: x and y coordinates of every intersection (columns 1 and 2), its ID (column 3) and whether it is a node on border (column 4), which means that this intersection is located at domain's boundary through which vehicles may enter (inflows), or exit (outflows);

2) "../ModelValidation/RoadTable.csv" -- contains information about roads: ID1 and ID2 (columns 3 and 4) are the id's of corresponding intersections that the road is connecting, ID_road (column 5) is the road's ID, max_vel (column 6) is its free-flow limit estimated from real measurements, then we have number of lanes (column 7) and road's length (column 8);

3) "../ModelValidation/RoadFRC.csv" -- contains information about road importance: column 1 is the road's ID and column 2 is the importance class from 1 to 7 with roads of class 1 being the most important.

4) "../ModelValidation/TurnTable.csv" -- contains turning ratios between any pair of roads: ID1 of incoming road (column 1), ID2 of outgoing road (column 2) and the turning ratio between these roads (column 5). 

II. Data from real sensors (model validation part):

4) "../ModelValidation/Timestamp.csv" -- contains time in seconds at which the data are given (unix timestamp), the time step equals to one minute;

5) "../ModelValidation/Density.csv" -- contains estimated density from real sensors: first number is road_id followed by its density (that is assumed to be constant within one road) at all time instants, then the next road_id with its density data for each time instant and so on;

6) "../ModelValidation/AllInflows.csv" -- contains inflow values (in veh/hour) for every road for every time step (one minute). If road is outgoing from intersection that is not on border, then the inflow value is zero;

7) "../ModelValidation/AllOutflows.csv" -- contains outflow values (in veh/hour) for every road for every time step (one minute). If road is incoming into intersection that is not on border, then the outflow value is zero.

The main file of the project is mainwindow.cpp: in its constructor we specify the file names to be loaded, start simulation starting time (line 26) and simulation step size (line 28). We can also change there the weighting parameter eta used to approximate parameters for every cell (line 4), and parameter d0 (line 5) is used for Gaussian Kernel estimation.

Other important classes are:
-- UrbanNetwork, which contains all the network geometry information (this is where all the network files are read). This network is used for both densities.
-- NSWEModel, which contains translation procedure of all network and intersection parameters into NSWE-formulation (function processIntersections). After all parameters are defined in NSWE, it calls constructInterpolation function that approximates these parameters defined for every intersection to be defined on every cell of a network. Then update is performed, where the Godunov numerical scheme is applied for the state update using NSWE model.
There is also a function getSSIMDiff_mean_weighted used to compute the weighted SSIM index between two densities (one from NSWE and the other one from real data).
-- GrenobleData, where all the data estimated from the real-life experiments are loaded. In function reconstructDensity the density initially given for each road is defined for every cell. Thereby, every road is divided in 10 parts and density values are presented as points on the border between these parts. Then Gaussian Kernel estimation is used to determine density for every cell in the domain.
-- TrafficSystem implements concurrent thread for parallel NSWE simulation relative to the main visualization thread