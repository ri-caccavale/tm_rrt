#include "../tm_rrt.h"

std::vector<PlanStepTM> TM_RRTplanner::plan_divided_BFS(State& start, State& goal, std::vector<PlanStepTM>& rrt_plan, double rrt_timeout, double horizon, double rrt_step, double low_level_rrt_timeout) {

    double stop_distance = 0.01;

    int n_samples = 0;
    int n_accomplished = 0;

    int n_samples_rejected = 0;

    bool plan_found = false;

    //std::vector<PlanStepTM> oldPlan = rrt_plan;
    //rrt_plan.clear();

    RRTree rrt;
    StepMT void_step;
    rrt.addNode(start, -1, 0, 100, void_step);

    std::cout << "INIT_state: " << std::endl;
    start.cout(false);
    std::cout << "GOAL_state: " << std::endl;
    goal.cout(false);

    RRTcluster rrt_cluster;
    std::vector<int> app = {0};
    rrt_cluster.createCluster(start.var, app);

    int par_id = 0;

    int first_to_delete = 0;



    std::cout << "done! " << std::endl;



    //INITIALIZATIONS - begin

    double coin;

    //SELECT RANDOM VAR-STATE
    double P_select_var = 0.5;
    double P_task_to_goal = 0.3;
    std::unordered_map < std::string, bool> v_random;

    //FIND NEAREST CLUSTER
    std::vector<int> vec_c_best;
    int best_c = -1;
    double best_c_dist = 1000;
    double c_dist = 0;

    //define the probability to select the goal as the new particle
    double P_go_to_goal = 0.3;

    //Task t_rand;
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
    

    std::cout << "RRT CLUSTER DIMENSIONS: " << rrt_cluster.size() << ", " << rrt_cluster.nodes[0].size() << ", " << rrt_cluster.vars[0].size() << std::endl;

    int not_exists_num = 0;
    int state_over_num = 0;

    

    //-- INITIALIZE task planner, tp_plans to reach the tp_states
    std::list< std::unordered_map < std::string, bool> > tp_states;
    std::list< std::list<Task> > tp_plans;

    std::list<Task> current_plan;

    double tp_starting_time = 0;
    double tp_elapsed_secs = 0;

    int tp_trial_counter = 0;

    std::cout<<"TP, randomizing tasks: "<<std::endl;
    std::vector<Task> original_task_set = task_set;
    task_set.clear();
    int r_task = 0;
    while(original_task_set.size()>0){
        r_task = iRand(0,original_task_set.size()-1);
        task_set.push_back(original_task_set[r_task]);
        original_task_set.erase(original_task_set.begin()+r_task);
    }
    std::cout<<"TP, plotting tasks: "<<std::endl;
    for(auto i=0; i<task_set.size(); i++)
        std::cout<<task_set[i].name<<std::endl;


    std::cout<<"TP, loading applicable tasks: "<<std::endl;
    //initialize open nodes
    for(auto i=0; i<task_set.size(); i++){  //normal order
    //for(int i=task_set.size()-1; i>=0; i--){  //inverted oreder
        if(is_applicable(task_set[i],start.var)){
            std::cout<<task_set[i].name<<" ";
            tp_states.push_back(start.var);
            std::list<Task> new_plan = {task_set[i]};
            tp_plans.push_back(new_plan);
        }
    }
    std::cout<<std::endl;
    //--

    
    if(debug_on)
        rrt_timeout = -1; //infinity

    //while you have time
    while ( ( rrt_timeout < 0 || elapsed_secs < rrt_timeout ) && !plan_found) {

        //-- BFS task planner
        bool tp_found = false; //find next solution

        //you should start again from skretch
        rrt_cluster.clear();
        rrt_cluster.createCluster(start.var, app);

        tp_starting_time = elapsed_secs;

        while( ( rrt_timeout < 0 || elapsed_secs < rrt_timeout ) && !tp_found && tp_states.size() > 0){

            timer.tic();

            std::unordered_map<std::string, bool> new_v_state = apply_to_state(tp_plans.front().back(), tp_states.front() );
            current_plan = tp_plans.front();

            tp_states.pop_front();
            tp_plans.pop_front();

            //for (auto it = current_plan.begin(); it != current_plan.end(); ++it){
            //    std::cout << it->name << " ";
            //}
            //std::cout<<std::endl;
            //sleep(1);

            if(distance(new_v_state, goal.var) < 0.1){
                tp_found  = true;

                tp_elapsed_secs += (elapsed_secs - tp_starting_time);

                std::cout<<"BFS-task-planner, plan found in "<<elapsed_secs - tp_starting_time<<" seconds"<<std::endl;

                //uncomment this to consider planning time also in the total time
                //elapsed_secs -= (elapsed_secs - tp_starting_time);

                continue;
                //now current_plan is the provided plan
            }

            for(auto i=0; i<task_set.size(); i++){
                if(is_applicable(task_set[i],new_v_state )){
                    tp_states.push_back(new_v_state);
                    //temporarly add the task
                    current_plan.push_back(task_set[i]);
                    tp_plans.push_back(current_plan);
                    //remove the task
                    current_plan.pop_back();

                }
            }
            //update the elapsed time
            elapsed_secs += timer.toc(); 
            
        }
        //--

        if(tp_found){
            for (auto it = current_plan.begin(); it != current_plan.end(); ++it){
                std::cout << it->name << " ";
            }
            std::cout<<std::endl;
            tp_trial_counter++;
        }
        else{
            std::cout<<"BFS-task-planner, plan NOT found!"<<std::endl;
            break;
        }



        //restart from initial state
        std::unordered_map<std::string, bool> v_current = start.var;
        int current_cluster = 0;

        //double low_level_rrt_timeout = 5.0; //1.0; //seconds
        double low_level_elapsed_secs = elapsed_secs + low_level_rrt_timeout;

        //-- RRT
        while( elapsed_secs < low_level_elapsed_secs && ( rrt_timeout < 0 || elapsed_secs < rrt_timeout ) && !plan_found) {

            timer.tic();

            n_samples++;

            //get task from the plan
            Task t_current = current_plan.front();

            //std::cout<<"rrt, current_task: "<<t_current.name<<std::endl;


            //SELECT RANDOM POSE

            //flip a float-coint between 0 and 1
            coin = fRand(0, 1);
            //if goal particle is selected
            if (coin <= P_go_to_goal){
                //set the random pose as the goal pose
                t_rand_goal = find_pose(v_current, t_current.target);
                //p_random = goal.pose;
                p_random = t_rand_goal;
            }
            //oth.
            else {
                p_random = Pose3d(fRand(-horizon, horizon), fRand(-horizon, horizon), fRand(-180, 180));
            }

            //create the random state
            S_random = State(v_current, p_random);

            // std::cout<<"RANDOM_state: "<<std::endl;
            // S_random.cout();

            seed::time::Clock tester;
            tester.tic();

            //std::cout<<"rrt, sampled "<<std::endl;

            //FIND NEAREST NODE OF THE RRT (here we use the weighted-distance!)
            S_near_index = 0;
            S_near_cluster = current_cluster; //0;
            S_near_best_distance = 10000;
            S_near_current_distance = -1;



            //now.. for each node of the cluster
            for (auto j = 0; j < rrt_cluster.nodes[current_cluster].size(); j++) {
                //compute the current distance (with the given symm_distance of the corresponding cluster)
                S_near_current_distance = distance(rrt.node[ rrt_cluster.nodes[current_cluster][j] ], S_random, 0.0);

                //if it is the first time, or the current distance is lower than the best one
                if (S_near_current_distance < S_near_best_distance) {
                    //select the current node of the cluster the nearest node
                    S_near = rrt.node[ rrt_cluster.nodes[current_cluster][j] ];
                    //remember, "nodes" is the vector of the indexes of the rrt
                    S_near_index = rrt_cluster.nodes[current_cluster][j];
                    //save also the cluster.. to speedup the insertion
                    S_near_cluster = current_cluster;
                    //update the best distance
                    S_near_best_distance = S_near_current_distance;
                }
            }

            //std::cout<<"rrt, node found "<<std::endl;
            

            if (S_near_best_distance >= 10000) {
                state_over_num++;
                continue; //the task is over
            }

            rrt_cluster.sampled[S_near_cluster]++;

            record_time += tester.toc();


            //NEW steer -> perform a path in the direction of p_random
            Pose3d task_target_pose = find_pose(S_near.var, t_current.target);
            std::vector<PlanStepTM> motion_plan = path_in_direction(S_near, S_random, task_target_pose, rrt_step); //works with 0.5m
            //steer - end
            
            //std::cout<<"rrt, steering done "<<std::endl;

            //if the path is consistent
            //if( motion_plan.size() != 0 ){
            if (motion_plan.size() > 0) { //non-trivial trajectory

                //std::cout<<"rrt, steering not trivial "<<std::endl;
                
                State S_final;

                for (int i = 0; i < motion_plan.size(); i++) {

                    StepMT step_new;

                    step_new.task = t_current;

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
                        n_accomplished++;
                        //add another node with the updated var state!
                        State S_new2 = S_new;
                        S_new2.var = apply_to_state(step_new.task, S_new);
                        
                        S_final = S_new2;

                        step_new.motion.fs = 0;
                        step_new.motion.ls = 0;
                        step_new.motion.ts = 0;

                        //double S_new2_cost = distance(S_new,S_near) + rrt.cost[S_near_index]; //DOES IT UPDATE COST?

                        gate_state_count[step_new.task.name]++;
                        
                        int new_id2 = rrt.addNode(S_new2, new_id, S_new_cost, S_new_obs, step_new);

                        //flag this task as accomplished in the current cluster
                        rrt_cluster.accomplished[S_near_cluster][step_new.task.name] = true;

                        if (distance(S_new2.var, goal.var) < 0.1) {

                            solution_time = elapsed_secs; //timer.toc();

                            std::cout << ansi::green << "RRT: solution found after " << solution_time << "secs (" << n_samples << " samples)" << ansi::end << std::endl;
                            //std::cout<<ansi::green<<"     solution obst: "<<S_new_obs<<" meters"<<ansi::end<<std::endl;
                            std::cout << ansi::green << "     solution cost: " << S_new_cost << " meters" << ansi::end << std::endl;
                            std::cout << ansi::green << "     solution rrt-id: " << new_id2 << " steps" << ansi::end << std::endl;
                            //cout_less_cluster(rrt_cluster,goal);

                            plan_found = true; //exit when the first plan is found!
                        }

                        //to plot
                        int old_size = rrt_cluster.size();
                        int new_cluster_index;

                        //add the new node to the associated cluster (NOTE: a new cluster can be created)
                        new_cluster_index = rrt_cluster.addToCluster(S_new2.var, new_id2);

                        if (old_size == new_cluster_index) {
                            if (distance(S_new2.var, goal.var) == 5)
                                std::cout << ansi::yellow << "NEW CLUSTER FOUND from " << S_near_cluster << " to " << new_cluster_index << " via " << step_new.task.name << "\t" << elapsed_secs << "\t" << distance(S_new2.var, goal.var) << ansi::end << std::endl;
                            else
                                std::cout << "NEW CLUSTER FOUND from " << S_near_cluster << " to " << new_cluster_index << " via " << step_new.task.name << "\t" << elapsed_secs << "\t" << distance(S_new2.var, goal.var) << std::endl;
                        }

                        //update the current symbolic state and remove the current task
                        v_current = S_new2.var;
                        
                        //rrt_cluster.cout(current_cluster);
                        current_cluster = new_cluster_index;

                        //std::cout<<t_current.name<<" done"<<std::endl;
                        current_plan.pop_front();

                        //std::cout<<"\t checking new action: "<< current_plan.front().name << " from State: "<< current_cluster <<std::endl;

                        //reset the RRT timeout
                        low_level_elapsed_secs = elapsed_secs + low_level_rrt_timeout;

                        break; //the rest of the path can be discard!
                    }                    
                    //add the new node to the best cluster (the var_state is the same by construction)
                    else if (i == motion_plan.size() - 1) //this looks faster but absolutely worst!!
                        rrt_cluster.addToCluster(S_near_cluster, new_id);

                }
                
                if(debug_on){
                    draw_tree(rrt,cv::Scalar(220,220,220));
                    std::vector<Point3d> void_obs;
                    draw_points(compute_dynamic_obstacles(S_near,void_obs), cv::Scalar(0,100,255) );
                    draw_sample(S_near,S_final,S_random);

                    wait_enter(ansi::green + "sample " + std::to_string(n_samples) + " accepted" + ansi::end);
                }
                
            } else {
                n_samples_rejected++;

                rrt_cluster.rejected_due_path[S_near_cluster]++;
                
                if(debug_on){
                    draw_tree(rrt,cv::Scalar(220,220,220));
                    std::vector<Point3d> void_obs;
                    draw_points(compute_dynamic_obstacles(S_near,void_obs), cv::Scalar(0,100,255) );
                    draw_sample(S_near,S_random);

                    wait_enter(ansi::red + "sample " + std::to_string(n_samples) + " rejected" + ansi::end);
                }
            }

            //update the elapsed time
            elapsed_secs += timer.toc();
        }
        // end RRT
    }

    if (!plan_found)
        std::cout << "NO SOLUTION FOUND IN " << elapsed_secs << " secs" << std::endl;
/*
    if(tp_trial_counter > 1){
        for(auto k=0; k<rrt_cluster.size(); k++){
            draw_map();
            draw_cluster_poses(rrt,rrt_cluster,k);
            rviz_image_plan();
            sleep(1);
            ros::spinOnce(); //probably useless
            char enter;
            std::cout<<"PRESS ENTER TO CONTINUE..."<<std::endl;
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
        }
    }
*/

    //this have to be done in two steps:
    //  find the best task then find the best state associated to the best task

    State S_final; //if we are lucky: p_final == goal
    double best_goal_dist = 10000;
    double best_task_dist = 10000;
    int index_final = 0;
    std::unordered_map < std::string, bool> task_var;

    //std::cout<<"DEFAULT BEST TASK: "<<rrt.act[0].task.name<<", task_dist: "<<best_task_dist<<", goal_dist: "<<best_goal_dist<<std::endl;

    //find the best task
    for (auto j = 0; j < rrt.size(); j++) {
        //NOTE: if the var-state is not reached, the pose distance have to be computed from the task-goal-pose!!
        double curr_task_dist, curr_goal_dist;

        if (rrt.act[j].task.name != "") {
            //task_var = apply_to_state( rrt.act[j].task, rrt.node[j] );
            //curr_task_dist = distance( task_var, goal.var );
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
            //std::cout<<"NEW BEST TASK FROM NODE "<<j<<" IS: "<<rrt.act[j].task.name<<", task_dist: "<<best_task_dist<<", goal_dist: "<<best_goal_dist<<std::endl;
        }
    }

    //std::cout<<"final index: "<<index_final<<std::endl;

    //find the path following the selected branch
    int k = index_final;

    std::vector<PlanStepTM *> plan_steps;

    double path_dist = 0; //to plot
    //for each parents until the start node
    while (k != -1) {
        //add the node to the plan
        int parent_index = rrt.parent[k];
        PlanStepTM *ps = new PlanStepTM();
        //ps.act = action_in_direction( rrt.node[parent_index], rrt.node[k] , goal);
        ps->act = rrt.act[k];
        //ps->state = rrt.node[parent_index];
        //ps->cost = rrt.cost[parent_index];
        //ps->min_obst = rrt.min_obst[parent_index];
        ps->state = rrt.node[k];
        ps->cost = rrt.cost[k];
        ps->min_obst = rrt.min_obst[k];

        if (parent_index >= 0)
            path_dist += distance2d(rrt.node[k].pose, rrt.node[parent_index].pose);

        //if(rrt.act[k].task.name == "")
        //    std::cout<<ansi::red<<"TASK IN THE PLAN IS NULL!!!"<<ansi::end<<std::endl;

        //rrt_plan.insert(rrt_plan.begin(), *ps);
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
        }
    }
/*
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
    //rrt_cluster.cout();
    //cout_cluster(rrt_cluster, goal);
    cout_cluster(rrt_cluster, goal, true); //only truthset
    std::cout << ansi::end;

    for (auto i = 0; i < rrt_cluster.size(); i++)
        if (distance(rrt_cluster.vars[i], goal.var) == 5) {
            std::cout << "PLOTTING: " << std::endl;
            rrt_cluster.cout(i);
            draw_cluster_poses(rrt, rrt_cluster, i);
            break;
        }
*/

    std::stringstream ss;
    ss << solution_time << "\t" << n_samples << "\t" << n_samples_rejected << "\t" << rrt_plan.size() << "\t" << path_dist << "\t" << rrt.cost[index_final] << "\t" << n_accomplished_in_plan << "\t" << tp_elapsed_secs; //<<"\t"<< record_time;
    rrt_report.push_back(ss.str());

    //return the plan
    return rrt_plan;

}
