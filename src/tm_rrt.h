#include <iostream>
#include <pthread.h>
#include <string>
#include <map>
#include <math.h>
#include <time.h>
#include <fstream>

#include "ros/ros.h"
#include "std_msgs/String.h"
#include "std_msgs/Float32.h"
#include "std_msgs/UInt32.h"

//#include "opencv/cv.h"
#include "opencv2/opencv.hpp"
#include "opencv2/ml/ml.hpp"
#include <random>
#include <execinfo.h>
#include <signal.h>
#include <ros/package.h>


#include "sensor_msgs/Joy.h"
#include "sensor_msgs/LaserScan.h"

#include "utils/seed_geometry.h"
#include "utils/seed_debug.h"
#include "utils/seed_time.h"

#include <string>

#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

#include <image_transport/image_transport.h>
#include <cv_bridge/cv_bridge.h>

#include <visualization_msgs/Marker.h>
#include <tf/transform_listener.h>
#include <tf/transform_broadcaster.h>

//#include "eclipseclp_interface/Eclipseclp_interface.h"
#include "swipl_interface/swipl_interface.h"

//#include "tm_rrt_vrep_sim.h"
//#include "tm_rrt_gazebo_sim.h"
//#include "tm_rrt_ompl.h"
//#include "tm_rrt_fcl.h"

//#include <time.h>

#define MOTIONPLANNER_STEP_FS 0.1 //m/s
#define MOTIONPLANNER_STEP_LS 0.1 //m/s
#define MOTIONPLANNER_STEP_TS 10 //deg/s
#define MOTIONPLANNER_PERIOD 0.5 //s

//sizes of robot
#define ROBOT_DEFAULT_LENGTH 0.8805
#define ROBOT_DEFAULT_WIDTH 0.55

#define ROBOT_CARRYING_LENGTH 1.25
#define ROBOT_CARRYING_WIDTH 0.94


//tmRRT constants
#define TM_RRT_WEIGHT_TOP 10
#define TM_RRT_WEIGHT_BOTTOM 1
#define TM_RRT_EXECUTING 0

#define RRT_DRAW_MAP_SIZE 1200//600

enum class TaskSelector{ BEST, UNIFORM, MONTECARLO };

//task represented in the planning domain
// pre-conditions/post_conditions are vectors of variables (eg. a -a) to be true

class Task {
public:
    //both are in {-1,0,1}
    Task();

    Task(std::string nname, std::vector<std::string>& npre, std::vector<std::string>& npost, std::string ntarget);

    // Overload = operator to set from Step3d.
    Task operator=(Task a);

    std::string name;
    std::vector<std::string> pre_conditions;
    std::vector<std::string> post_conditions;
    std::string target;
    Pose3d pose;

};

class Step3d {
public:
    //both are in {-1,0,1}
    Step3d();

    Step3d(double nfs, double nls, deg180 nts);

    // Overload = operator to set from Step3d.
    Step3d operator=(Step3d a);

    double fs; //m/s
    double ls; //m/s
    deg180 ts; //deg/s
};

class StepMT {
public:
    //both are in {-1,0,1}
    StepMT();

    StepMT(double nfs, double nls, deg180 nts, Task& ntask);

    StepMT(Step3d nmotion, Task& ntask);

    // Overload = operator to set from Step3d.

    inline StepMT operator=(StepMT a) {
        this->motion = a.motion;
        this->task = a.task;
        return a;
    }

    void cout();

    Step3d motion;
    Task task;
};

class State {
public:
    State();

    State(std::unordered_map<std::string, bool>& nvariables, Pose3d npose);

    inline State operator=(State a) {
        this->var = a.var;
        this->pose = a.pose;

        return a;
    }

    void cout(bool true_only);

    std::string toString(bool true_only, bool symbolic_only);

    //std::vector<bool> var;
    std::unordered_map<std::string, bool> var;
    Pose3d pose;

};

class PlanStepTM {
public:
    PlanStepTM();

    PlanStepTM(State& nstate, StepMT& nact);

    PlanStepTM(State& nstate, Step3d nmotion, Task& ntask);

    PlanStepTM(Pose3d npos, std::unordered_map<std::string, bool>& nvar_set, Step3d nmotion, Task& ntask);

    PlanStepTM(State& nstate, StepMT& nact, double ncost, double nobst);

    PlanStepTM(Pose3d npos, std::unordered_map<std::string, bool>& nvar_set, Step3d nmotion, Task& ntask, double ncost, double nobst);

    inline PlanStepTM operator=(PlanStepTM a) {
        this->state = a.state;

        this->act = a.act;

        this->cost = a.cost;
        this->min_obst = a.min_obst;
        return a;
    }

    void cout();

    //state
    State state;
    //step
    StepMT act;
    double cost; //cost from the start point
    double min_obst; //cummulative distance from obstacles
};

class RRTree {
public:
    RRTree();

    RRTree(int n);

    RRTree(RRTree& nrrt);

    RRTree operator=(RRTree a);

    int addNode(State& nnode, int nid, double ncost, double nobs, StepMT& nact);

    inline int size() {
        return node.size();
    }

    //get the path starting from the slected index (node)
    std::vector<PlanStepTM> path(int index);

    //all vectors share the same index
    std::vector<State> node; //explored nodes
    std::vector<double> cost; //costs from the start point
    std::vector<int> parent; //parents of the nodes
    std::vector<double> min_obst; //total obstacle distance
    std::vector<StepMT> act; //step which lead to the cuurent state from the previous one
};

//this groups together the nodes with the same var-state
//  NOTE: cluster is misleading ...better situation? frame? circumstance? eventuality? 

class RRTcluster {
public:
    RRTcluster();

    inline int createCluster(std::unordered_map<std::string, bool>& nvar, std::vector<int>& nnode){
        vars.push_back( nvar );
        nodes.push_back( nnode );
        
        accomplished.push_back( std::unordered_map<std::string,bool>() );
        
        sampled.push_back(0);
        rejected_due_task.push_back(0);
        rejected_due_path.push_back(0);
        
        number_of_clusters++;

        //return nodes.size()-1;
        return number_of_clusters-1;
    }

    inline int addToCluster(std::unordered_map<std::string, bool>& nvar, int nid, bool autoCreate = true) {

        int c = findCluster(nvar);

        if (c >= 0)
            nodes[c].push_back(nid);
        else if (autoCreate) {
            std::vector<int> app(1, nid);
            c = createCluster(nvar, app);
        }

        return c;
    }

    inline int addToCluster(int c, int nid) {

        if (c <= nodes.size())
            nodes[c].push_back(nid);
        else
            return -1;

        return c;
    }
    int findCluster(std::unordered_map<std::string, bool>& nvar);

    inline int size(){
        //return nodes.size();
        return number_of_clusters;
    }

    inline void clear(){
        vars.clear();
        nodes.clear();
        accomplished.clear();
        sampled.clear();
        rejected_due_task.clear();
        rejected_due_path.clear();
        number_of_clusters = 0;
    }

    void cout();
    
    void cout(int _id);

    //plots just the id and the number of nodes
    void cout_less();

    //both share the same index
    //  associate the explored var_state to the explored RRT nodes (by index)
    std::vector< std::unordered_map<std::string, bool> > vars;
    std::vector< std::vector<int> > nodes;
    //the initial state is the first explored
    
    //optimization: skip accomplished tasks! -> not working in general, it depends on the starting node of the forest!! 
    std::vector< std::unordered_map<std::string, bool> > accomplished; //list of the accomplished tasks of the cluster
    
    //debug, number of times this cluster is sampled
    std::vector< int > sampled;
    std::vector< int > rejected_due_task;
    std::vector< int > rejected_due_path;
    
    int number_of_clusters;
};

class TM_RRTplanner {
public:
    
    TM_RRTplanner(std::string path_to_node_directory, std::string domain_file_name) ;

    //plan = plan_rrt_simple(S_init,S_goal,plan,60,5,path_len); //example of a call

    void set_goal_state(std::unordered_map<std::string, bool> &tS_goal, Pose3d &bS_goal);

    void set_initial_state(std::unordered_map<std::string, bool> &tS_init, Pose3d &bS_init);

    void set_initial_state(Pose3d &bS_init);

    //rrt_timeout is secs
    //std::vector<PlanStepTM> plan_rrt_simple(double rrt_timeout = 0.1, double horizon = 2, double rrt_step = 0.5);
    std::vector<PlanStepTM> plan_rrt_simple(State& start, State& goal, std::vector<PlanStepTM>& rrt_plan, double rrt_timeout = 0.1, double horizon = 2, double rrt_step = 0.5);

    std::vector<PlanStepTM> plan_rrt_naive(State& start, State& goal, std::vector<PlanStepTM>& rrt_plan, double rrt_timeout = 0.1, double horizon = 2, double rrt_step = 0.5);

    std::vector<PlanStepTM> plan_divided_BFS(State& start, State& goal, std::vector<PlanStepTM>& rrt_plan, double rrt_timeout = 0.1, double horizon = 2, double rrt_step = 0.5, double low_level_rrt_timeout = 5.0);


    std::vector< PlanStepTM > transform_plan(std::vector< PlanStepTM > vec, std::string in_frame, std::string out_frame);


    void rviz_plot_plan(std::vector< PlanStepTM > vec);

    void rviz_image_plan();
    
    void draw_map();

    void draw_map(Pose3d pose, Pose3d goal);

    void draw_points(Pose3d pose, std::vector<Point3d> p, cv::Scalar col = cv::Scalar(255, 255, 255));
    
    void draw_points(std::vector<Point3d> p, cv::Scalar col = cv::Scalar(255, 255, 255), double m2px = 100, double offset = RRT_DRAW_MAP_SIZE / 2);

    void draw_plan(Pose3d pose, int start_index = 1);

    void draw_plan(std::vector<PlanStepTM> in_plan, Pose3d pose);

    void draw_robot_pose(double px, double py, double pw, double r_len, double r_wid, cv::Scalar color, double m2px = 100, double offset = RRT_DRAW_MAP_SIZE / 2);

    void draw_robot_pose(Pose3d p, double r_len, double r_wid, cv::Scalar color, double m2px = 100, double offset = RRT_DRAW_MAP_SIZE / 2);

    void draw_robot_pose(State state, cv::Scalar color, double m2px = 100, double offset = RRT_DRAW_MAP_SIZE / 2);

    void draw_robot_pose(double px, double py, double pw, State state, cv::Scalar color, double m2px = 100, double offset = RRT_DRAW_MAP_SIZE / 2);

    void draw_path(std::vector<PlanStepTM> p, cv::Scalar color, double m2px = 100, double offset = RRT_DRAW_MAP_SIZE / 2);

    void draw_state(State s, cv::Scalar color, double m2px = 100, double offset = RRT_DRAW_MAP_SIZE / 2);
    
    void draw_sample(State s_near, State s_new, State s_sample);
    
    void draw_sample(State s_near, State s_sample);
    
    void draw_tree(RRTree rrt, cv::Scalar color, double m2px = 100, double offset = RRT_DRAW_MAP_SIZE / 2);
    
    // on(c2,p3) draw_cluster_poses(RRTree t, RRTcluster c, "on(c2,p3)")
    //returns the id of the cluster
    int draw_cluster_poses(RRTree t, RRTcluster c, std::vector<std::string> var);

    void draw_cluster_poses(RRTree t, RRTcluster c, int id);

    void cout_cluster(RRTcluster c, State _goal, bool truthset_only = false);

    //lpots only id, number of nodes and distance
    void cout_less_cluster(RRTcluster c, State _goal);

//private:
    
    Pose3d find_pose(std::unordered_map<std::string, bool>& var, std::string toFind);

    std::vector<Task> tasks_from_PROLOG(std::vector<std::string>& task_name);

    void domain_from_PROLOG();

    bool is_pose_reached(Pose3d p1, Pose3d p2);

    std::string get_object_loacation(State &state, std::string object_to_find);

    std::vector<Point3d> compute_dynamic_obstacles(State &state);

    std::vector<Point3d> compute_dynamic_obstacles(State &state, std::vector<Point3d> &obs);

    //WARNING: this function is context specific!!!
    std::vector<Point3d> add_object_to_obstacles(std::string obj_name, State &state);

    //WARNING: this function is context specific!!!
    std::vector<Point3d> add_object_to_obstacles(std::string obj_name, State &state, std::vector<Point3d> &obs);

    std::vector<Point3d> remove_carrying_from_obstacles(State &state, std::vector<Point3d> &obs, double consis_margin = 0.02);

    //WARNING: this function is context specific!!!
    std::vector<Point3d> remove_object_from_obstacles(std::string obj_name, State &state, std::vector<Point3d> &obs, double consis_margin = 0.02);

    /** check if the pose is consistent adapting robot len considering the state
     *
     * @param state: the state to be checked pose and variables
     * @return 0 if the pose is colliding, otherwise the distance with the closer
     *      obstacle is returned.
     * 
     * //WARNING: this function is context specific!!!
     */
    double is_pose_consistent(State &state, double consis_margin = 0.02);

    /** check if the pose is consistent adapting robot len considering the state
     *
     * @param state: the state to be checked pose and variables
     * @param obstacles: the obstacles you can find in the current state
     * @return 0 if the pose is colliding, otherwise the distance with the closer
     *      obstacle is returned.
     * 
     * //WARNING: this function is context specific!!!
     */
    double is_pose_consistent(State &state, std::vector<Point3d> &updated_obstacles, double consis_margin = 0.02);

    /** check if the pose is consistent (ie. not colliding obstacles)
     *
     * @param pose: the pose to be checked (x,y,yaw)
     * @param r_len, r_wid: robot length and robot width respectively
     * @return 0 if the pose is colliding, otherwise the distance with the closer
     *      obstacle is returned.
     */
    double is_pose_consistent(Pose3d pose, double r_len, double r_wid, std::vector<Point3d> &obs, double consis_margin = 0.02);

    //distance between poses (3d)

    inline double distance(Pose3d p1, Pose3d p2) {
        //Euclidean distance for the linear components (x,y)
        double lin_dist = sqrt(((p1.x - p2.x)*(p1.x - p2.x))+((p1.y - p2.y)*(p1.y - p2.y)));

        //arc of the robot rotation
        double ang_dist = dtor(std::abs(p1.w - p2.w)) * robot_length / 2;

        if (dtor(std::abs(p1.w - p2.w)) > 180)
            std::cout << "DISTANCE-ERROR ang dist was " << dtor(std::abs(p1.w - p2.w)) << std::endl;

        //sum both distances
        return lin_dist + ang_dist;
    }

    //distance between states-variable-vector (var, motion)
    //  NOTE: the maps are passed as reference! this avoid the copy which is time consuming

    inline double distance(std::unordered_map<std::string, bool>& v1, std::unordered_map<std::string, bool>& v2) {
        //complexity is n*m

        if (v1.size() > v2.size())
            return symmetric_difference(v2, v1);
        else
            return symmetric_difference(v1, v2);
    }


    double symmetric_difference(std::unordered_map<std::string, bool>& v1, std::unordered_map<std::string, bool>& v2);

    //distance between states (var, motion)

    inline double distance(State& s1, State& s2) {

        //double w_b = 1;//1; //bottom -> for motion
        //double w_t = 10; //top    -> for variables

        //weighted-sum [ISSUE: when the goal-pose is close, it can win over the task!]
        double d = w_b * distance(s1.pose, s2.pose) + w_t * distance(s1.var, s2.var);

        //power
        //double d = distance(s1.pose,s2.pose) * std::pow( 10 , distance(s1.var, s2.var) );

        //sum both distances
        return d;
    }

    //distance between states (optimized)
    //  NOTE: you can provide the symmetric distance to avoid computation (time optimization)

    inline double distance(State& s1, State& s2, double symm_distance) {

        //double w_b = 1;//1; //bottom -> for motion
        //double w_t = 10; //top    -> for variables

        //weighted-sum [ISSUE: when the goal-pose is close, it can win over the task!]
        double d = w_b * distance(s1.pose, s2.pose) + w_t * symm_distance;

        //power
        //double d = distance(s1.pose,s2.pose) * std::pow( 10 , distance(s1.var, s2.var) );

        //sum both distances
        return d;
    }

    //simple euclidean distance (skip the yaw)

    inline double distance2d(Pose3d p1, Pose3d p2) {
        return sqrt(((p1.x - p2.x)*(p1.x - p2.x))+((p1.y - p2.y)*(p1.y - p2.y)));
    }

    //compute the line between 2 poses, returns a, b and c of the line

    inline void getLine(Pose3d l1, Pose3d l2, double &a, double &b, double &c) {
        a = l1.y - l2.y;
        b = l2.x - l1.x;
        c = l1.x * l2.y - l2.x * l1.y;
    }

    //compute the distance from p to the line touching l1 and l2

    inline double distance_from_line(Pose3d l1, Pose3d l2, Pose3d p) {
        //compute the line
        double a, b, c;
        getLine(l1, l2, a, b, c);
        //return std::abs(a * p.x + b * p.y + c) / std::sqrt(a * a + b * b);
        //if the 2 points are the same
        if (a * a + b * b == 0)
            return 0;
        //otherwise return the distance
        return (a * p.x + b * p.y + c) / std::sqrt(a * a + b * b);
    }

    void update_obstacles_from_LIDAR();

    void update_obstacles_from_FILE(std::string path_to_map = "");

    //transform a vector of points from robot to world frame
    std::vector<Point3d> robot2world_frame(std::vector<Point3d> vec);

    //OPTIMIZATION of RRT: find shortcuts bypassing noisy paths
    std::vector<PlanStepTM> random_shortcut_plan(std::vector<PlanStepTM>& in_plan, double rs_timeout = 0.01);

    //get the plan-step that is closer to the current robot pose
    int get_nearest_plan_step();

    bool check_plan();

    std::unordered_map<std::string, bool> apply_to_state(Task& t, State& s);

    std::unordered_map<std::string, bool> apply_to_state(Task& t, std::unordered_map<std::string, bool>& v);

    bool is_applicable(Task& t, State& s);

    bool is_applicable(Task& t, std::unordered_map<std::string, bool>& inv);

    StepMT action_in_direction(State& s1, State& s2, State& g);

    Task task_in_direction(State& s1, State& s2, TaskSelector mode = TaskSelector::BEST);

    Task task_in_direction(std::unordered_map<std::string, bool>& v1, std::unordered_map<std::string, bool>& v2, TaskSelector mode = TaskSelector::BEST);
    
    Task task_in_direction_BEST(State& s1, State& s2);

    Task task_in_direction_BEST(std::unordered_map<std::string, bool>& v1, std::unordered_map<std::string, bool>& v2);

    Task task_in_direction_UNIF(State& s1, State& s2);

    Task task_in_direction_UNIF(std::unordered_map<std::string, bool>& v1, std::unordered_map<std::string, bool>& v2);

    Task task_in_direction_MC(State& s1, State& s2);

    Task task_in_direction_MC(std::unordered_map<std::string, bool>& v1, std::unordered_map<std::string, bool>& v2);

    
    Step3d step_in_direction(State& s1, State& s2, State& g);

    //steer functions (RRT)
    Step3d step_in_direction(Pose3d p1, Pose3d p2, Pose3d goal);

    //states have to be the same vars!!
    std::vector< PlanStepTM > path_in_direction(State& s_start, State& s_stop, Pose3d& p_goal, double stop_distance = 1);

    Pose3d estimate_new_pose(Pose3d p, Step3d s);

    inline double fRand(double fMin, double fMax) {
        // double f = (double)rand() / RAND_MAX;
        // return fMin + f * (fMax - fMin);

        std::uniform_real_distribution<double> unif(fMin, fMax);
        // std::default_random_engine re;
        return unif(random_generator); //unif(re);
    }

    inline int iRand(int iMin, int iMax) {
        // return rand() % iMax + iMin;
        std::uniform_int_distribution<int> unif(iMin, iMax);
        // std::default_random_engine re;
        return unif(random_generator); //unif(re);
    }

    double median(std::vector<double> scores);

    inline void wait_enter(std::string to_plot){
        std::string enter;
        std::cout<<to_plot<<"...";
        std::getline(std::cin,enter);
    }

    
//attributes

    
    swipl_interface *swi;

    //the "w" com is 0
    //std::vector<Pose3d> obstacles;
    std::vector<Point3d> obstacles;
    std::vector<Point3d> laser_points;

    //partial path generated by low-leve planner (mainly used for debug)
    std::vector<PlanStepTM> partial_path; 

    //constants
    double robot_length;
    double robot_width;

    double x;
    double y;
    deg180 yaw;

    //variables considered during planning
    std::vector<std::string> var_set;
    //map of poses
    std::unordered_map<std::string, Pose3d> pose_map;
    //tasks considered during planning
    std::vector<Task> task_set;

    //goal state
    State S_goal;
    //initial state
    State S_init;

    std::vector<PlanStepTM> plan;

    int curr_plan_step;
    bool going_to_goal;

    //DEBUG
    bool once; //for debug
    bool debug_on; //for debug
    bool cout_debug_on; //for debug
    std::vector<std::string> rrt_report;
    int path_looping_count;
    double record_time;

    cv::Mat image;

    //TM-RRT PARAMETERS
    double w_b; //weight of the task (top)
    double w_t; //weight of the path (bottom)
    double path_len; //max length of the path

    
    //std::default_random_engine re;
    std::mt19937 random_generator; //Standard mersenne_twister_engine seeded with rd()
    
    
    //OPTIMIZATIONS
    bool second_chance_heuristic;

    
    ros::NodeHandle nh;
    ros::Publisher marker_pub;
    image_transport::Publisher image_pub;
    

    tf::TransformListener laser_listener;
    std::vector< PlanStepTM > plan_odo_frame;
    

    tf::TransformListener target_listener;
    tf::TransformBroadcaster target_broadcaster;
    

    tf::StampedTransform transform;

    //simulation
    std::unordered_map<std::string, Object2d> simulation;

    std::string map_file;
    std::string path_to_directory;
    
    TaskSelector mode;
};