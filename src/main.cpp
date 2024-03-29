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
    
    //get prameters from ROS

    //DOMAIN FILE IS NECESSARY
    std::string param_domain_file;
    if( !ros::param::get("/tm_rrt/domain_file", param_domain_file) ){
        std::cout<<ansi::red<<"ERROR: planning domain not specified"<<ansi::end<<std::endl;
        return 0;
    }

    //get debug value 
    bool param_debug = get_ros_param<bool>("/tm_rrt/debug",false);
    //get planner type (simple, divided)
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
    //go-to-goal probabilities
    double param_p_s = get_ros_param<double>("/tm_rrt/p_s",0.3);
    double param_p_c = get_ros_param<double>("/tm_rrt/p_c",0.3);

    //load symbolic doman and geometric map 

    TM_RRTplanner tmRRT(ros::package::getPath("tm_rrt"), "domains/" + param_domain_file);
    
    tmRRT.update_obstacles_from_FILE();

    //ROS-based drawing and visualization of the map
    
    tmRRT.draw_map(tmRRT.S_init.pose,tmRRT.S_goal.pose);
        
    std::vector<Point3d> void_obs;
    tmRRT.draw_points(tmRRT.S_init.pose, tmRRT.compute_dynamic_obstacles(tmRRT.S_init,void_obs) );
    
    tmRRT.rviz_image_plan();

    int j = 1;

    //execute it several times
    for (auto i = 0; i < param_number_of_runs; i++) {

        tmRRT.plan.clear();
        tmRRT.draw_map();

        std::cout << std::endl << "TM_PLANNER: ROUND NUMBER " << i + 1 << std::endl << std::endl;

        //set parameters
        tmRRT.w_b = param_w_b; //weight of the path (bottom) 
        tmRRT.w_t = param_w_t; //weight of the task (top)
        tmRRT.path_len = param_path_len; //max length of the path (default value)

        tmRRT.P_task_to_goal = param_p_s; //probability to sample the goal symbolic state
        tmRRT.P_go_to_goal = param_p_c; //probability to sample the goal configuration

        tmRRT.debug_on = param_debug;
        
        tmRRT.mode = TaskSelector(param_task_selector);

        if(param_planner_type == "tm_rrt") //TM-RRT
            tmRRT.plan_TM_RRT(tmRRT.S_init, tmRRT.S_goal, tmRRT.plan, param_planner_timeout, param_rrt_horizon, tmRRT.path_len);
        else if(param_planner_type == "bfs_rrt") //BFS+RRT
            tmRRT.plan_BFS_RRT(tmRRT.S_init, tmRRT.S_goal, tmRRT.plan, param_planner_timeout, param_rrt_horizon, tmRRT.path_len, param_bsf_timeout);
        else{
            std::cout<<ansi::red<<"ERROR: planner of type "<<param_planner_type<<" does not exists"<<ansi::end<<std::endl;
            return 0;
        }


        cv::putText(tmRRT.image, std::to_string(j) + "." + std::to_string(i + 1), cv::Point(15, 30), cv::FONT_HERSHEY_DUPLEX, 1, cv::Scalar(0, 143, 143), 2);
        tmRRT.draw_plan(tmRRT.S_init.pose);
        tmRRT.rviz_image_plan();
        
        if(!ros::ok())
            break;

        tmRRT.ros_publish_plan();
        ros::spinOnce();
    }

    tmRRT.rrt_report.push_back("");

    std::cout << std::endl << "REPORT: " << std::endl;
    std::cout << "\t" << "planning_time" << "\t" << "explored_states" << "\t" << "rejected_states" << "\t" << "plan_size" << "\t" << "plan_length" << "\t" << "combined_cost" << "\t" << "number_of_tasks"<<"\t"<< "BFS_time"<<std::endl;
    for (auto i = 0; i < tmRRT.rrt_report.size(); i++)
        std::cout << i + 1 << "\t" << tmRRT.rrt_report[i] << std::endl;
        
    return 0;
}
