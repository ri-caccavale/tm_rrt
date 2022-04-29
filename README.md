# tm_rrt testing package
This package has been designed to assess the performance of the TM-RRT (Task and Motion planner based on Rapidly-exploring Random Trees) considering several problem configurations in hospital-logistic scenario with a mobile robotic platform.

The symbolic knowledge along with the problems formulations are available as Prolog files in the ```domains/``` folder

The source code of the planner, based on a naive version of the RRT, is available in the ```src/planners/planner_simple.cpp``` file

For comparison pourposes a divided planner based on simple BFS (symbolic) and RRT (geometric) is available in the ```src/planners/divided_BFS.cpp``` file

The code is wrapped in a ROS node allowing fast launch and parameter setting. The output of the node is a table in the terminal showing the performance of the selected planner. In addition a 2D graphical representation of generated plans is available from the /image_plan topic.

## Installation requirements
Install and configure ROS from [here](http://wiki.ros.org/ROS/Installation)

Install SWI-Prolog from repository:
```
sudo apt-add-repository ppa:swi-prolog/stable
sudo apt-get update
sudo apt-get install swi-prolog
```

## Running instructions
Run the test node using the following command:
```
roslaunch tm_rrt tm_rrt.launch
```

## Reference
Please acknowledge the authors in any acedemic publication that used parts of these codes.
```
@inproceedings{caccavale2021combining,
  title={Combining Task and Motion Planning through Rapidly-Exploring Random Trees},
  author={Caccavale, Riccardo and Finzi, Alberto},
  booktitle={2021 European Conference on Mobile Robots (ECMR)},
  pages={1--6},
  organization={IEEE}
}
```
