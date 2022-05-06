#include "../tm_rrt.h"

std::vector<PlanStepTM> TM_RRTplanner::plan_TM_RRT(State& start, State& goal, std::vector<PlanStepTM>& rrt_plan, double rrt_timeout, double horizon, double rrt_step) {

    double stop_distance = 0.01;

    int n_samples = 0;
    int n_accomplished = 0;

    int n_samples_rejected = 0;

    bool plan_found = false;

    RRTree rrt;
    StepMT void_step;

    //init the root node
    rrt.addNode(start, -1, 0, 100, void_step);

    std::cout << "INIT_state: " << std::endl;
    start.cout(false);
    std::cout << "GOAL_state: " << std::endl;
    goal.cout(false);

    //NOTE: clusters are used to group all RRT nodes having the same symbolic state
    //      this is a simple optimization to speed-up the search inside the trees
    RRTcluster rrt_cluster;
    std::vector<int> app = {0};
    rrt_cluster.createCluster(start.var, app);

    int par_id = 0;

    int first_to_delete = 0;

    std::cout << "done! " << std::endl;



    //INITIALIZATIONS - begin

    double coin;

    //SELECT RANDOM VAR-STATE
    double P_select_var = 0.5; //probability to select true or false symbolic variable
    std::unordered_map < std::string, bool> v_random;

    //FIND NEAREST CLUSTER
    std::vector<int> vec_c_best;
    int best_c = -1;
    double best_c_dist = 1000;
    double c_dist = 0;

    Task t_rand;
    Pose3d t_rand_goal;

    Pose3d p_random;
    State S_random;

    //FIND NEAREST NODE OF THE RRT (here we use the weighted-distance!)
    State S_near;
    int S_near_index = 0;
    int S_near_cluster = 0;
    double S_near_best_distance = -1;
    double S_near_current_distance = -1;
    double c_symm_dist = 0;

    //INITIALIZATIONS - end

    std::unordered_map<std::string, int> gate_state_count;

    //seed clock!
    record_time = 0;
    double elapsed_secs = 0;
    double solution_time = 0;
    seed::time::Clock timer;
    timer.tic();

    std::cout << "RRT CLUSTER DIMENSIONS: " << rrt_cluster.size() << ", " << rrt_cluster.nodes[0].size() << ", " << rrt_cluster.vars[0].size() << std::endl;

    int not_exists_num = 0;
    int state_over_num = 0;
    
    if(debug_on)
        rrt_timeout = -1; //infinity

    //while you have time
    while ( ( rrt_timeout < 0 || elapsed_secs < rrt_timeout ) && !plan_found) {
        n_samples++;

        //---------- RANDOM_STATE ----------//

        //SELECT RANDOM VAR-STATE (RANDOM_SYM)
        v_random.clear();

        coin = fRand(0, 1);
        if (coin <= P_task_to_goal) {
            v_random = goal.var;

            going_to_goal = true;
        } else {
            //for each possible variable
            for (auto i = 0; i < var_set.size(); i++) {
                //flip a coin
                coin = fRand(0, 1);
                //if good, select the variable for the rand state
                if (coin <= P_select_var) {
                    v_random[var_set[i]] = true;
                }
            }

            going_to_goal = false;
        }

        //FIND NEAREST CLUSTER (NEAREST_SYM)
        vec_c_best.clear();
        best_c = -1;
        best_c_dist = 1000;
        c_dist = 0;

        //find all min-distance clusters
        for (auto i = 0; i < rrt_cluster.size(); i++) {
            c_dist = distance(v_random, rrt_cluster.vars[i]);

            if (c_dist < best_c_dist) {
                best_c = i;
                best_c_dist = c_dist;

                vec_c_best.clear();
                vec_c_best.push_back(i);
            } else if (c_dist == best_c_dist) {
                vec_c_best.push_back(i);
            }
        }

        //select randomly the best cluster among those with the same (min) distance
        best_c = vec_c_best[ iRand(0, vec_c_best.size() - 1) ];

        //FIND THE ACTION (SELECT_ACTION)

        //find best task from the best cluster_var to the rand_var
        Task t_rand = task_in_direction(rrt_cluster.vars[best_c], v_random, mode); //selector

        //if the selected task not exists
        if (t_rand.name == "") {
            not_exists_num++;
            continue;
        }

        //SELECT RANDOM CONFIGURATION AND GOAL (RAND_CONF)

        //flip a float-coint between 0 and 1
        coin = fRand(0, 1);
        //if goal particle is selected
        if (coin <= P_go_to_goal){
            //set the random pose as the goal pose (apply G function)
            t_rand_goal = find_pose(rrt_cluster.vars[best_c], t_rand.target);
            p_random = t_rand_goal;
        }
        //oth.
        else {
            p_random = Pose3d(fRand(-horizon, horizon), fRand(-horizon, horizon), fRand(-180, 180));
        }

        //create the random combined state
        S_random = State(v_random, p_random);

        if (t_rand.name == "") {
            State toPlot(rrt_cluster.vars[best_c], p_random);
            std::cout << ansi::red << "NO task available for state: " << ansi::end << std::endl;
            toPlot.cout(false);
        }

        // std::cout<<"RANDOM_state: "<<std::endl;
        // S_random.cout();

        //---------- EXTEND ----------//

        seed::time::Clock tester;
        tester.tic();

        //FIND NEAREST NODE OF THE RRT (here we use the weighted-distance!) (NEAREST_NEIGHBOR)
        S_near_index = 0;
        S_near_cluster = 0;
        S_near_best_distance = 10000;
        S_near_current_distance = -1;
        //for each node
        // NOTE: we will explore the rrt cluster-by-cluster so to compute the symmetric distance (complex) only once for each cluster!
        c_symm_dist = 0;
        //so... for each cluster
        for (auto i = 0; i < rrt_cluster.size(); i++) {

            //compute symm_distance between the random state and the variables in the cluster
            c_symm_dist = distance(S_random.var, rrt_cluster.vars[i]);

            //branch-and-bound -like simple optimization
            if(c_symm_dist > S_near_best_distance)
                continue;

            //skip tasks that are not applicable in this cluster (reject less samples)
            if (!is_applicable(t_rand, rrt_cluster.vars[i]))
                continue;

            //now.. for each node of the cluster
            for (auto j = 0; j < rrt_cluster.nodes[i].size(); j++) {
                //compute the current distance (with the given symm_distance of the corresponding cluster)
                S_near_current_distance = distance(rrt.node[ rrt_cluster.nodes[i][j] ], S_random, c_symm_dist);

                //if it is the first time, or the current distance is lower than the best one
                if (S_near_current_distance < S_near_best_distance) {
                    //select the current node of the cluster the nearest node
                    S_near = rrt.node[ rrt_cluster.nodes[i][j] ];
                    //remember, "nodes" is the vector of the indexes of the rrt
                    S_near_index = rrt_cluster.nodes[i][j];
                    //save also the cluster.. to simplify the insertion
                    S_near_cluster = i;
                    //update the best distance
                    S_near_best_distance = S_near_current_distance;
                }
            }
        }

        //in case no best distance is found
        if (S_near_best_distance >= 10000) {
            state_over_num++;
            continue; //the task is over
        }

        //counter for statistics
        rrt_cluster.sampled[S_near_cluster]++;

        record_time += tester.toc();

        //FIND MOTION in the random direction (NEW_STATE)
        //      NOTE: here we use the previously sampled symbolic action 

        //STEERING -> perform a path in the direction of p_random
        Pose3d task_target_pose = find_pose(S_near.var, t_rand.target);
        std::vector<PlanStepTM> motion_plan = path_in_direction(S_near, S_random, task_target_pose, rrt_step);

        //if the path is consistent (CHECK)
        if (motion_plan.size() > 0) { //non-trivial trajectory
            
            State S_final;

            //add the segment to the RRT (ADD_VERTEX + ADD_EDGE)

            for (int i = 0; i < motion_plan.size(); i++) {

                StepMT step_new;

                step_new.task = t_rand;

                step_new.motion = motion_plan[i].act.motion;

                double obs_dist = motion_plan[i].min_obst;

                State S_new;
                S_new.var = S_near.var;

                S_new.pose = motion_plan[i].state.pose; //= estimate_new_pose(S_near.pose,step_new.motion); //compute new position

                S_final = S_new;
                
                int new_id;
                //compute the lenght of the RRT branch (meters)
                double S_new_cost = distance(S_new, S_near) + rrt.cost[S_near_index];
                //update the closest obstacle of the branch until now
                double S_new_obs = obs_dist < rrt.min_obst[S_near_index] ? obs_dist : rrt.min_obst[S_near_index];

                // //push back p_new directly
                new_id = rrt.addNode(S_new, S_near_index, S_new_cost, S_new_obs, step_new);

                S_near_index = new_id;

                //CHECK if subgoal is reached
                if (is_pose_reached(S_new.pose, find_pose(S_new.var, step_new.task.target))) {
                    // here we set the new symbolic state to the final step of the trajectory
                    n_accomplished++;
                    //add another node with the updated var state!
                    State S_new2 = S_new;
                    S_new2.var = apply_to_state(step_new.task, S_new);
                    
                    S_final = S_new2;

                    step_new.motion.fs = 0;
                    step_new.motion.ls = 0;
                    step_new.motion.ts = 0;

                    gate_state_count[step_new.task.name]++;
                    
                    int new_id2 = rrt.addNode(S_new2, new_id, S_new_cost, S_new_obs, step_new);

                    //flag this task as accomplished in the current cluster (for statistics)
                    rrt_cluster.accomplished[S_near_cluster][step_new.task.name] = true;

                    if (distance(S_new2.var, goal.var) < 0.1) {

                        solution_time = timer.toc();

                        std::cout << ansi::green << "RRT: solution found after " << solution_time << "secs (" << n_samples << " samples)" << ansi::end << std::endl;
                        std::cout << ansi::green << "     solution cost: " << S_new_cost << " meters" << ansi::end << std::endl;
                        std::cout << ansi::green << "     solution rrt-id: " << new_id2 << " steps" << ansi::end << std::endl;

                        plan_found = true; //exit when the first plan is found!
                    }

                    //to plot
                    int old_size = rrt_cluster.size();
                    int new_cluster_index;

                    //add the new node to the associated cluster (NOTE: a new cluster can be created)
                    new_cluster_index = rrt_cluster.addToCluster(S_new2.var, new_id2);

                    if (old_size == new_cluster_index) {
                        if (distance(S_new2.var, goal.var) == 5)
                            std::cout << ansi::yellow << "NEW CLUSTER FOUND from " << S_near_cluster << " to " << new_cluster_index << " via " << step_new.task.name << "\t" << timer.toc() << "\t" << distance(S_new2.var, goal.var) << ansi::end << std::endl;
                        else
                            std::cout << "NEW CLUSTER FOUND from " << S_near_cluster << " to " << new_cluster_index << " via " << step_new.task.name << "\t" << timer.toc() << "\t" << distance(S_new2.var, goal.var) << std::endl;
                    }

                    break;
                }                    
                //add the new node to the best cluster (the var_state is the same by construction)
                else if (i == motion_plan.size() - 1)
                    rrt_cluster.addToCluster(S_near_cluster, new_id);
            }
            
            if(debug_on){
                draw_map();
                draw_path(partial_path, cv::Scalar(0,100,0));
                draw_cluster_poses(rrt,rrt_cluster,S_near_cluster);
                std::vector<Point3d> void_obs;
                draw_points(compute_dynamic_obstacles(S_near,void_obs), cv::Scalar(0,100,255) );
                draw_sample(S_near,S_final,S_random);
                rviz_image_plan();
                ros::spinOnce();
                wait_enter(ansi::green + "sample " + std::to_string(n_samples) + " accepted" + ansi::end);
            }
            
        } else {
            n_samples_rejected++;

            rrt_cluster.rejected_due_path[S_near_cluster]++;
            
            if(debug_on){
                draw_map();
                draw_path(partial_path, cv::Scalar(0,0,100));
                draw_cluster_poses(rrt,rrt_cluster,S_near_cluster);
                std::vector<Point3d> void_obs;
                draw_points(compute_dynamic_obstacles(S_near,void_obs), cv::Scalar(0,100,255) );
                draw_sample(S_near,S_random);
                rviz_image_plan();
                ros::spinOnce();
                wait_enter(ansi::red + "sample " + std::to_string(n_samples) + " rejected" + ansi::end);
            }
        }

        //update the elapsed time
        elapsed_secs = timer.toc();

    }

    // FIND SOLUTION: EXTRACT THE PALN FROM THE RRT 

    if (!plan_found)
        std::cout << "NO SOLUTION FOUND IN " << elapsed_secs << " secs" << std::endl;

    //this have to be done in two steps:
    //  find the best task then find the best state associated to the best task

    State S_final;
    double best_goal_dist = 10000;
    double best_task_dist = 10000;
    int index_final = 0;
    std::unordered_map < std::string, bool> task_var;

    //find the best task
    for (auto j = 0; j < rrt.size(); j++) {
        double curr_task_dist, curr_goal_dist;

        if (rrt.act[j].task.name != "") {
            curr_task_dist = distance(rrt.node[j].var, goal.var);
            curr_goal_dist = distance(rrt.node[j].pose, find_pose(rrt.node[j].var, rrt.act[j].task.target));
        } else
            continue;

        if ((curr_task_dist == best_task_dist && curr_goal_dist < best_goal_dist) ||
                (curr_task_dist < best_task_dist)) {

            S_final = rrt.node[j];
            index_final = j;
            best_goal_dist = curr_goal_dist;
            best_task_dist = curr_task_dist;
        }
    }

    //find the path following the selected branch
    int k = index_final;

    std::vector<PlanStepTM *> plan_steps;

    double path_dist = 0;
    //for each parents until the start node
    while (k != -1) {
        //add the node to the plan
        int parent_index = rrt.parent[k];
        PlanStepTM *ps = new PlanStepTM();
        ps->act = rrt.act[k];
        ps->state = rrt.node[k];
        ps->cost = rrt.cost[k];
        ps->min_obst = rrt.min_obst[k];

        if (parent_index >= 0)
            path_dist += distance2d(rrt.node[k].pose, rrt.node[parent_index].pose);

        plan_steps.push_back(ps);

        k = parent_index;
    }

    int n_accomplished_in_plan = 0;
    int cd = 0;
    std::cout << "RRT-plan:" << std::endl;
    for (auto i = 0; i < plan_steps.size(); i++) {
        cd = (plan_steps.size() - 1) - i;
        rrt_plan.push_back(*plan_steps[cd]);
        if (plan_steps[cd]->act.motion.fs == 0 && plan_steps[cd]->act.motion.ls == 0 && plan_steps[cd]->act.motion.ts == 0){
            n_accomplished_in_plan++;
            std::cout << "\t " << plan_steps[cd]->act.task.name << std::endl;
            std::cout << "\t\t "; plan_steps[cd]->state.cout(true);
        }
    }

    std::cout << ansi::cyan;
    std::cout << "RRT-report" << std::endl;
    std::cout << "\t final index: " << index_final << std::endl;
    std::cout << "\t samples: " << n_samples << std::endl;
    for(auto it = gate_state_count.begin(); it != gate_state_count.end(); it++)
        std::cout << "\t\t gates of " << it->first << ": "<< it->second << std::endl;
    std::cout << "\t rejected: " << n_samples_rejected << std::endl;
    std::cout << "\t rejected (not-exists): " << not_exists_num << std::endl;
    std::cout << "\t rejected (state-over): " << state_over_num << std::endl;
    std::cout << "\t looping paths: " << path_looping_count << std::endl;
    std::cout << "\t rrt-size: " << rrt.size() << std::endl;
    std::cout << "\t plan-size: " << rrt_plan.size() << std::endl;
    std::cout << "\t subtasks: " << n_accomplished << std::endl;
    std::cout << "\t solution cost: " << rrt.cost[index_final] << std::endl;
    std::cout << "\t solution euclid dist: " << path_dist << std::endl;
    std::cout << "\t solution dist-to-goal: " << distance(rrt.node[index_final], goal) << std::endl;
    std::cout << "\t solution euclid dist-to-goal: " << distance(rrt.node[index_final].pose, goal.pose) << std::endl;
    std::cout << "EXPLORED clusters:" << std::endl;
    std::cout << ansi::end;

    std::stringstream ss;
    ss << solution_time << "\t" << n_samples << "\t" << n_samples_rejected << "\t" << rrt_plan.size() << "\t" << path_dist << "\t" << rrt.cost[index_final] << "\t" << n_accomplished_in_plan; //<<"\t"<< record_time;
    rrt_report.push_back(ss.str());

    for (auto i = 0; i < rrt_cluster.size(); i++)
        if (distance(rrt_cluster.vars[i], goal.var) == 5) {
            std::cout << "PLOTTING: " << std::endl;
            rrt_cluster.cout(i);
            draw_cluster_poses(rrt, rrt_cluster, i);
            break;
        }

    return rrt_plan;
}
