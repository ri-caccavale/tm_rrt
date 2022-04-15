#include "tm_rrt.h"

template <typename T> 
T get_ros_param(std::string param_name, T default_value){
    T param_value;
    if( !ros::param::param<T>(param_name, param_value, default_value) )
        std::cout<<ansi::yellow<<"WARNING: parameter "<<param_name<<" not set, default value \""<<default_value<<"\" will be used"<<ansi::end<<std::endl;
    return param_value;
}

int main(int argc, char** argv) {
    
    ros::init(argc, argv, "tm_rrt");
    
//    //VrepSimulation vrep;
//    GazeboSimulation gazebo;
//    gazebo.test(argc, argv);
//    std::cout<<"after sim creation"<<std::endl;
//    return 0;
    
    //FCLsimulation fcls;
    //fcls.test<double>(argc, argv);
    //std::cout<<"after sim creation"<<std::endl;
    //return 0;
    
    //TM_RRTplanner tmRRT(ros::package::getPath("tm_rrt"), "PLAN_nofree_1c2p.prolog"); //case1
    //TM_RRTplanner tmRRT(ros::package::getPath("tm_rrt"), "PLAN_nofree_2c4p.prolog"); //case2
    //TM_RRTplanner tmRRT(ros::package::getPath("tm_rrt"), "PLAN_nofree_3c6p.prolog"); //case3
    //TM_RRTplanner tmRRT(ros::package::getPath("tm_rrt"), "PLAN_nofree_4c8p.prolog"); //case4

    //TM_RRTplanner tmRRT(ros::package::getPath("tm_rrt"), "PLAN_nofree_2c4p (anomaly2.0).prolog"); //case2.0
    //TM_RRTplanner tmRRT(ros::package::getPath("tm_rrt"), "PLAN_nofree_2c4p (anomaly2.1).prolog"); //case2.1
    //TM_RRTplanner tmRRT(ros::package::getPath("tm_rrt"), "PLAN_nofree_3c6p (anomaly2.2).prolog"); //case2.2
    //TM_RRTplanner tmRRT(ros::package::getPath("tm_rrt"), "PLAN_nofree_3c6p (anomaly2.25).prolog"); //case2.25
    //TM_RRTplanner tmRRT(ros::package::getPath("tm_rrt"), "PLAN_nofree_3c6p (anomaly2.3).prolog"); //case2.3
    //TM_RRTplanner tmRRT(ros::package::getPath("tm_rrt"), "PLAN_nofree_4c8p (anomaly2.4).prolog"); //case2.4 //TOO HARD
    //TM_RRTplanner tmRRT(ros::package::getPath("tm_rrt"), "PLAN_nofree_4c8p (anomaly2.5).prolog"); //case2.5

    /* OLD
    //TM_RRTplanner tmRRT(ros::package::getPath("tm_rrt"), "PLAN_nofree_2c5p (anomaly2.1).prolog"); //case2.1
    //TM_RRTplanner tmRRT(ros::package::getPath("tm_rrt"), "PLAN_nofree_3c6p (anomaly2.2).prolog"); //case2.2
    //TM_RRTplanner tmRRT(ros::package::getPath("tm_rrt"), "PLAN_nofree_3c6p (anomaly2.3).prolog"); //case2.3
    //TM_RRTplanner tmRRT(ros::package::getPath("tm_rrt"), "PLAN_nofree_3c6p (anomaly2.4).prolog"); //case2.4
    */
    
    //TM_RRTplanner tmRRT(ros::package::getPath("tm_rrt"), "PLAN_nofree_2c3p (PlanRob).prolog"); //case5

    //TM_RRTplanner tmRRT(ros::package::getPath("tm_rrt"), "PLAN_nofree_2c4p (local minima).prolog"); //case6
    //TM_RRTplanner tmRRT(ros::package::getPath("tm_rrt"), "PLAN_nofree_3c6p (anomaly2).prolog"); //case7

    //TM_RRTplanner tmRRT(ros::package::getPath("tm_rrt"), "PLAN_nofree_4c7p (anomaly3).prolog"); //case8

    //get prameters from ROS

    //DOMAIN FILE IS NECESSARY
    std::string param_domain_file;
    if( !ros::param::get("/tm_rrt/domain_file", param_domain_file) ){
        std::cout<<ansi::red<<"ERROR: planning domain not specified"<<ansi::end<<std::endl;
        return 0;
    }

    //get debug value 
    bool param_debug = get_ros_param<bool>("/tm_rrt/debug",false);
    //get planner type (simple, naive, divided)
    std::string param_planner_type = get_ros_param<std::string>("/tm_rrt/planner_type","simple");
    //get task selector (best = 0, uniform = 1, montecarlo = 2)
    int param_task_selector = get_ros_param<int>("/tm_rrt/task_selector",0);
    //get number of runs for this domain
    int param_number_of_runs = get_ros_param<int>("/tm_rrt/number_of_runs",30);
    //get timeout for the planner (in seconds)
    double param_planner_timeout = get_ros_param<double>("/tm_rrt/planner_timeout",300);
    //get timeout for the BFS (in seconds)
    double param_bsf_timeout = get_ros_param<double>("/tm_rrt/BFS_timeout",2);
    //get horizon of the rrt (in meters)
    double param_rrt_horizon = get_ros_param<double>("/tm_rrt/RRT_horizon",5);
    //get TM-RRT parameters
    double param_w_b = get_ros_param<double>("/tm_rrt/w_b",1.0);
    double param_w_t = get_ros_param<double>("/tm_rrt/w_t",5.0);
    double param_path_len = get_ros_param<double>("/tm_rrt/path_len",0.9);


    TM_RRTplanner tmRRT(ros::package::getPath("tm_rrt"), "domains/" + param_domain_file);

    
    tmRRT.update_obstacles_from_FILE();
    
    tmRRT.draw_map(tmRRT.S_init.pose,tmRRT.S_goal.pose);
        
    std::vector<Point3d> void_obs;
    tmRRT.draw_points(tmRRT.S_init.pose, tmRRT.compute_dynamic_obstacles(tmRRT.S_init,void_obs) );
    
    tmRRT.rviz_image_plan();
       
    //            for(auto j=1; j<=10; j++){
    //            
    //            std::cout<<std::endl<<"TM_PLANNER: TRAIN NUMBER "<<j<<std::endl;
    //            //path_len = 0.5 + (0.1*j);
    //            //std::cout<<std::endl<<"\t path_len: "<<path_len<<std::endl<<std::endl;
    //            
    //            w_b = 1.0;
    //            w_t = w_b * (j*0.5);
    //            std::cout<<std::endl<<"\t w_b: "<<w_b<<" w_t: "<<w_t<<" (rate: "<<(j*0.5)<<")"<<std::endl<<std::endl;
    int j = 1;

    //execute it several times
    for (auto i = 0; i < 30; i++) {

        tmRRT.plan.clear();
        tmRRT.draw_map();

        std::cout << std::endl << "TM_PLANNER: ROUND NUMBER " << i + 1 << std::endl << std::endl;

        //tmRRT.second_chance_heuristic = true;

        //set parameters
        //tmRRT.w_b = 1.0; //weight of the path (bottom) 
        //tmRRT.w_t = 5.0; //weight of the task (top)
        //tmRRT.path_len = 0.9; //max length of the path (default value)
        tmRRT.w_b = param_w_b; //weight of the path (bottom) 
        tmRRT.w_t = param_w_t; //weight of the task (top)
        tmRRT.path_len = param_path_len; //max length of the path (default value)

        //tmRRT.debug_on = true;
        tmRRT.debug_on = param_debug;
        
        //tmRRT.mode = TaskSelector::BEST;
        //tmRRT.mode = TaskSelector::UNIFORM;
        //tmRRT.mode = TaskSelector::MONTECARLO;
        tmRRT.mode = TaskSelector(param_task_selector);

        //tmRRT.plan = tmRRT.plan_rrt_simple(tmRRT.S_init, tmRRT.S_goal, tmRRT.plan, 300, 5, tmRRT.path_len); //TM-RRT
        //tmRRT.plan = tmRRT.plan_rrt_naive(tmRRT.S_init, tmRRT.S_goal, tmRRT.plan, 300, 5, tmRRT.path_len);
        //tmRRT.plan = tmRRT.plan_divided_BFS(tmRRT.S_init, tmRRT.S_goal, tmRRT.plan, 300, 5, tmRRT.path_len, 2.0);
        if(param_planner_type == "simple")
            tmRRT.plan_rrt_simple(tmRRT.S_init, tmRRT.S_goal, tmRRT.plan, param_planner_timeout, param_rrt_horizon, tmRRT.path_len);
        else if(param_planner_type == "divided")
            tmRRT.plan_divided_BFS(tmRRT.S_init, tmRRT.S_goal, tmRRT.plan, param_planner_timeout, param_rrt_horizon, tmRRT.path_len, param_bsf_timeout);
        else{
            std::cout<<ansi::red<<"ERROR: planner of type "<<param_planner_type<<" does not exists"<<ansi::end<<std::endl;
            return 0;
        }


        cv::putText(tmRRT.image, std::to_string(j) + "." + std::to_string(i + 1), cv::Point(15, 30), cv::FONT_HERSHEY_DUPLEX, 1, cv::Scalar(0, 143, 143), 2);
        tmRRT.draw_plan(tmRRT.S_init.pose);
        tmRRT.rviz_image_plan();

        //std::cout << "obstacles size: " << tmRRT.obstacles.size() << std::endl;
        
        if(!ros::ok())
            break;
    }

    tmRRT.rrt_report.push_back("");

    std::cout << std::endl << "RRT REPORT: " << std::endl;
    for (auto i = 0; i < tmRRT.rrt_report.size(); i++)
        std::cout << i + 1 << "\t" << tmRRT.rrt_report[i] << std::endl;

    //            }


    //std::cout << "TM_PLANNER: optimize " << std::endl;
    //optimize
    //plan = random_shortcut_plan(plan);

    /*
    std::cout << "PLOT-PLAN: " << std::endl;
    for (auto i = 0; i < tmRRT.plan.size(); i++) {
        if (tmRRT.plan[i].act.task.name == "") {
            std::cout << ansi::red;
            tmRRT.plan[i].cout();
            std::cout << ansi::end;
        } else {
            std::cout << i << "- ";
            tmRRT.plan[i].cout();
        }

        std::cout << "->" << std::endl;
    }
    */

        
    return 0;
}
