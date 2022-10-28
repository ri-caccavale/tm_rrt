# tm_rrt testing package
This package has been designed to assess the performance of the TM-RRT (Task and Motion planner based on Rapidly-exploring Random Trees) considering several problem configurations in hospital-logistic scenario with a mobile robotic platform.

The symbolic knowledge along with an example of problem formulation is available as Prolog file in the ```domains/``` folder.

The source code of the planner, based on a naive version of the RRT, is available in the ```src/planners/planner_TM_RRT.cpp``` file.

For comparison pourposes a 2-layer planner based on simple BFS (symbolic) and RRT (geometric) is available in the ```src/planners/planner_BFS_RRT.cpp``` file.

The code is wrapped into a ROS node allowing fast launch and parameter setting. The output of the node is a table in the terminal showing the performance of the selected planner. In addition a 2D graphical representation of generated plans is available from the /image_plan topic while the list of the actions is published on the /tm_rrt_plan topic.

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
@article{caccavale2022rapidly,
  title={A rapidly-exploring random trees approach to combined task and motion planning},
  author={Caccavale, Riccardo and Finzi, Alberto},
  journal={Robotics and Autonomous Systems},
  volume={157},
  pages={104238},
  year={2022},
  publisher={Elsevier}
}
```
