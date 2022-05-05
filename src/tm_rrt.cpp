
#include "tm_rrt.h"

//    ***** Task *****    //

Task::Task() {
    name = "";
}

Task::Task(std::string nname, std::vector<std::string>& npre, std::vector<std::string>& npost, std::string ntarget) {
    this->name = nname;
    this->pre_conditions = npre;
    this->post_conditions = npost;
    this->target = ntarget;
}
// Overload = operator to set from Step3d.

Task Task::operator=(Task a) {
    this->name = a.name;
    this->pre_conditions = a.pre_conditions;
    this->post_conditions = a.post_conditions;
    this->target = a.target;

    return a;
}


//    ***** Step3d *****    //

Step3d::Step3d() {
    fs = 0;
    ls = 0;
    ts = 0;
};

Step3d::Step3d(double nfs, double nls, deg180 nts) {
    fs = nfs * MOTIONPLANNER_STEP_FS;
    ls = nls * MOTIONPLANNER_STEP_LS;
    ts = nts * MOTIONPLANNER_STEP_TS;
};

// Overload = operator to set from Step3d.

Step3d Step3d::operator=(Step3d a) {
    this->fs = a.fs;
    this->ls = a.ls;
    this->ts = a.ts;
    return a;
}


//    ***** StepMT *****    //

StepMT::StepMT() {
    motion.fs = 0;
    motion.ls = 0;
    motion.ts = 0;
}

StepMT::StepMT(double nfs, double nls, deg180 nts, Task& ntask) {
    motion.fs = nfs * MOTIONPLANNER_STEP_FS;
    motion.ls = nls * MOTIONPLANNER_STEP_LS;
    motion.ts = nts * MOTIONPLANNER_STEP_TS;

    task = ntask;
}

StepMT::StepMT(Step3d nmotion, Task& ntask) {
    motion = nmotion;
    task = ntask;
}

void StepMT::cout() {
    std::cout << "ACT( " << task.name << " " << motion.fs << " " << motion.ls << " " << motion.ts << " )." << std::endl;
}


//    ***** State *****    //

State::State() {
}

State::State(std::unordered_map<std::string, bool>& nvariables, Pose3d npose) {
    this->var = nvariables;
    this->pose = npose;
}

void State::cout(bool true_only) {
    std::cout << "STATE( ( ";
    for (auto it : var) {
        if(true_only && it.second)
            std::cout << it.first << " ";
        else if(!true_only)
            std::cout << it.first << ":" << it.second << " ";
    }
    std::cout << "), ( " << pose.x << " " << pose.y << " " << pose.w << " ) )." << std::endl;
}

std::string State::toString(bool true_only, bool symbolic_only){
    
    std::stringstream ss;
    
    bool first = true;
    
    ss << "state( ( ";
    for (auto it : var) {
        
        if(!first)
            ss << ", ";
        
        if(true_only && it.second){
            ss << it.first ;
            first = false;
        }
        else if(!true_only){
            first = false;
            
            if(it.second)
                ss << it.first ;
            else
                ss << "-" << it.first ;
        }
    }
    
    if(symbolic_only)
        ss << ") )";
    else
        ss << "), ( " << pose.x << " " << pose.y << " " << pose.w << " ) )";
    
    return ss.str();
}


//    ***** PlanStepTM *****    //

PlanStepTM::PlanStepTM() {
    state = State();

    act = StepMT();

    cost = 0;
    min_obst = 0;
}

PlanStepTM::PlanStepTM(State& nstate, StepMT& nact) {
    state = nstate;

    act = nact;

    cost = 0;
    min_obst = 0;
}

PlanStepTM::PlanStepTM(State& nstate, Step3d nmotion, Task& ntask) {
    state = nstate;

    act.motion = nmotion;
    act.task = ntask;

    cost = 0;
    min_obst = 0;
}

PlanStepTM::PlanStepTM(Pose3d npos, std::unordered_map<std::string, bool>& nvar_set, Step3d nmotion, Task& ntask) {
    state.pose = npos;
    state.var = nvar_set;

    act.motion = nmotion;
    act.task = ntask;

    cost = 0;
    min_obst = 0;
}

PlanStepTM::PlanStepTM(State& nstate, StepMT& nact, double ncost, double nobst) {
    state = nstate;

    act = nact;

    cost = ncost;
    min_obst = nobst;
}

PlanStepTM::PlanStepTM(Pose3d npos, std::unordered_map<std::string, bool>& nvar_set, Step3d nmotion, Task& ntask, double ncost, double nobst) {
    state.pose = npos;
    state.var = nvar_set;

    act.motion = nmotion;
    act.task = ntask;

    cost = ncost;
    min_obst = nobst;
}

void PlanStepTM::cout() {
    std::cout << "state( ( ";
    for (auto it : state.var) {
        std::cout << it.first << ":" << it.second << " ";
    }
    std::cout << "), ( " << state.pose.x << " " << state.pose.y << " " << state.pose.w << " ) )." << std::endl;

    std::cout << "act( " << act.task.name << " " << act.motion.fs << " " << act.motion.ls << " " << act.motion.ts << " )." << std::endl;
}

std::string PlanStepTM::toString() {
    //std::cout << "state( ( ";
    //for (auto it : state.var) {
    //    std::cout << it.first << ":" << it.second << " ";
    //}
    //std::cout << "), ( " << state.pose.x << " " << state.pose.y << " " << state.pose.w << " ) )." << std::endl;
    std::stringstream ss;
    ss << "( " << act.task.name << ", (" << act.motion.fs << ", " << act.motion.ls << ", " << act.motion.ts << " ) )";
    return ss.str();
}


//    ***** RRTree *****    //

RRTree::RRTree() {
    node = std::vector<State>();
    cost = std::vector<double>();
    parent = std::vector<int>();
    min_obst = std::vector<double>();
    act = std::vector<StepMT>();
}

RRTree::RRTree(int n) {
    node = std::vector<State>(n);
    cost = std::vector<double>(n);
    parent = std::vector<int>(n);
    min_obst = std::vector<double>(n);
    act = std::vector<StepMT>(n);
}

RRTree::RRTree(RRTree& nrrt) {
    node = nrrt.node;
    cost = nrrt.cost;
    parent = nrrt.parent;
    min_obst = nrrt.min_obst;
    act = nrrt.act;
}

RRTree RRTree::operator=(RRTree a) {
    node = a.node;
    cost = a.cost;
    parent = a.parent;
    min_obst = a.min_obst;
    act = a.act;

    return a;
}

int RRTree::addNode(State& nnode, int nid, double ncost, double nobs, StepMT& nact) {

    node.push_back(nnode);
    cost.push_back(ncost);
    min_obst.push_back(nobs);
    parent.push_back(nid);
    act.push_back(nact);

    return node.size() - 1;
}

//get the path starting from the slected index (node)

std::vector<PlanStepTM> RRTree::path(int index) {

    std::vector<PlanStepTM> rrt_path;

    //find the path following the selected branch
    int k = index;
    //add the best final node to the path
    PlanStepTM ps;
    ps.state = node[k];
    ps.cost = cost[k];
    ps.min_obst = min_obst[k];
    rrt_path.insert(rrt_path.begin(), ps);
    //for each parents until the start node
    while (k != 0) {
        //add the node to the plan
        int parent_index = parent[k];
        PlanStepTM ps;
        //ps.act = action_in_direction( rrt.node[parent_index], rrt.node[k] , goal);
        ps.act = act[k];
        ps.state = node[parent_index];
        ps.cost = cost[parent_index];
        ps.min_obst = min_obst[parent_index];

        rrt_path.insert(rrt_path.begin(), ps);
        k = parent_index;
    }

    return rrt_path;
}


//    ***** RRTcluster *****    //

RRTcluster::RRTcluster(){
    vars.clear();
    nodes.clear();

    accomplished.clear();

    sampled.clear();
    rejected_due_task.clear();
    rejected_due_path.clear();

    number_of_clusters = 0;
}

int RRTcluster::findCluster(std::unordered_map<std::string, bool>& nvar) {
    bool diff;
    //search for a cluster with the target var_map
    for (auto i = 0; i < vars.size(); i++) {
        //for each cluster, it is not different by dafault
        diff = false;
        //for each element of the cluster var_map
        for (auto it : vars[i]) {
            //if the element is not present in the in_var_map
            if (nvar[it.first] != it.second) {
                //the cluster is different
                diff = true;
                break;
            }
        }
        //if the cluster is not yet different
        if (!diff) {
            //for each element of the in_var_map
            for (auto it : nvar) {
                //if the element is not present in the cluster var_map
                if (vars[i][it.first] != it.second) {
                    //che cluster is different
                    diff = true;
                    break;
                }
            }
            //if the two maps are equivalent
            if (!diff) {
                //returns the cluster
                return i;
            }
        }
    }
    // this means that the cluster not exists
    return -1;
}

void RRTcluster::cout() {
    for (int i = 0; i < vars.size(); i++) {
        std::cout << "CLUSTER" << i << "( ( ";
        for (auto it : vars[i]) {
            std::cout << it.first << ":" << it.second << " ";
        }
        std::cout << ")" << std::endl;

        std::cout << "\t nodes in the cluster: " << nodes[i].size() << std::endl;
    }
}

void RRTcluster::cout(int _id){
    std::cout<<"CLUSTER"<<_id<<"( ( ";
    for(auto it : vars[_id]){
        std::cout<<it.first<<":"<<it.second<<" ";
    }
    std::cout<<")"<<std::endl;

    std::cout<<"\t nodes in the cluster: "<<nodes[_id].size()<<std::endl;
}


//plots just the id and the number of nodes
void RRTcluster::cout_less() {
    for (int i = 0; i < vars.size(); i++) {
        std::cout << "CLUSTER" << i << ": " << nodes[i].size() << std::endl;
    }
}


//    ***** TM_RRTplanner *****    //

TM_RRTplanner::TM_RRTplanner(std::string path_to_node_directory, std::string domain_file_name) : random_generator(std::random_device()()) {

    image_transport::ImageTransport it(nh);
    image_pub = it.advertise("image_plan", 1);

    plan_pub = nh.advertise<std_msgs::String>("tm_rrt_plan", 1);

    x = 0;
    y = 0;
    yaw = 0;

    once = true;

    cout_debug_on = false;
    debug_on = false;

    path_looping_count = 0;

    image = cv::Mat(1000, 1000, CV_8UC3, cv::Scalar(0, 0, 0));

    swi = new swipl_interface();
    swi->consult(path_to_node_directory + "/" + domain_file_name);
    
    path_to_directory = path_to_node_directory;
    
    curr_plan_step = 0;

    
    second_chance_heuristic = false;
    
    
    domain_from_PROLOG();

    std::cout << "TM_PLANNER: ...done " << std::endl;

    std::cout << std::endl << "TM_PLANNER: init done... " << std::endl;
    std::cout << "\t var_set size: " << var_set.size() << std::endl;
    std::cout << "\t task_set size: " << task_set.size() << std::endl;
}

void TM_RRTplanner::set_goal_state(std::unordered_map<std::string, bool> &tS_goal, Pose3d &bS_goal) {
    S_goal = State(tS_goal, bS_goal);
}

void TM_RRTplanner::set_initial_state(std::unordered_map<std::string, bool> &tS_init, Pose3d &bS_init) {
    S_init = State(tS_init, bS_init);
}

void TM_RRTplanner::set_initial_state(Pose3d &bS_init) {
    S_init.pose = bS_init;
}

// ************ SWI Prolog functions ************ //

std::vector<Task> TM_RRTplanner::tasks_from_PROLOG(std::vector<std::string>& task_name) {

    std::vector<Task> tsk(task_name.size());

    //per ogni schema del piano
    for (int i = 0; i < task_name.size(); i++) {

        //LOAD TASK NAME
        tsk[i].name = task_name[i];

        std::cout << "\t " << task_name[i] << " pre_cond:" << std::endl;
        //LOAD PRE_COND
        //tsk[i].pre_conditions = eclipse->query("getConstraint(" + (task_name[i]) + ")");
        tsk[i].pre_conditions = swipl_interface::list2vector(swipl_interface::functor2vector(swi->query("precondition(" + (task_name[i]) + ",X)"))[2]);

        std::cout << "\t " << task_name[i] << " post_cond:" << std::endl;
        //LOAD POST_COND
        //tsk[i].post_conditions = eclipse->query("getEffect(" + (task_name[i]) + ")");
        tsk[i].post_conditions = swipl_interface::list2vector(swipl_interface::functor2vector(swi->query("effect(" + (task_name[i]) + ",X)"))[2]);

        std::cout << "\t " << task_name[i] << " target:" << std::endl;
        //LOAD TARGET
        //tsk[i].target = eclipse->query("getTarget(" + (task_name[i]) + ")")[0];
        tsk[i].target = swipl_interface::list2vector(swipl_interface::functor2vector(swi->query("target(" + (task_name[i]) + ",X)"))[2])[0];

        std::cout << "\t ----------" << std::endl;
    }

    return tsk;

}

void TM_RRTplanner::domain_from_PROLOG() {

    //load VARIABLES:

    //var_set = eclipse->query("getVarList");
    var_set = swipl_interface::list2vector(swipl_interface::functor2vector(swi->query("findall(V,variable(V),L)"))[3]);

    //load TASKS:

    //std::vector<std::string> task_name = eclipse->query("getTaskList");
    std::vector<std::string> task_name = swipl_interface::list2vector(swipl_interface::functor2vector(swi->query("findall(T,operator(T),L)"))[3]);

    for (int i = 0; i < task_name.size(); i++) {

        Task tsk;

        //LOAD TASK NAME
        tsk.name = task_name[i];

        std::cout << "\t " << task_name[i] << " pre_cond:" << std::endl;
        //LOAD PRE_COND
        //tsk.pre_conditions = eclipse->query("getConstraint(" + (task_name[i]) + ")");
        tsk.pre_conditions = swipl_interface::list2vector(swipl_interface::functor2vector(swi->query("precondition(" + (task_name[i]) + ",X)"))[2]);

        std::cout << "\t " << task_name[i] << " post_cond:" << std::endl;
        //LOAD POST_COND
        //tsk.post_conditions = eclipse->query("getEffect(" + (task_name[i]) + ")");
        tsk.post_conditions = swipl_interface::list2vector(swipl_interface::functor2vector(swi->query("effect(" + (task_name[i]) + ",X)"))[2]);

        std::cout << "\t " << task_name[i] << " target:" << std::endl;
        //LOAD TARGET
        //tsk.target = eclipse->query("getTarget(" + (task_name[i]) + ")")[0];
        tsk.target = swipl_interface::list2vector(swipl_interface::functor2vector(swi->query("target(" + (task_name[i]) + ",X)"))[2])[0];

        std::cout << "\t ----------" << std::endl;

        task_set.push_back(tsk);
    }

    //load SHAPES of objects:

    //std::vector<std::string> objects = eclipse->query("getObjectList");
    std::vector<std::string> objects = swipl_interface::list2vector(swipl_interface::functor2vector(swi->query("findall(O,object(O,S),L)"))[3]);

    for (int i = 0; i < objects.size(); i++) {

        Object2d obj;

        //std::vector<std::string> shapes = eclipse->query("getObjectShapes(" + objects[i] + ")");
        std::vector<std::string> shapes = swipl_interface::list2vector(swipl_interface::functor2vector(swi->query("object(O,S)"))[2]);

        for (int j = 0; j < shapes.size(); j++) {

            std::vector<std::string> shape_vec = swipl_interface::functor2vector(shapes[j]);

            obj.addShapeRect(ston(shape_vec[1]),
                             ston(shape_vec[2]),
                             ston(shape_vec[3]),
                             ston(shape_vec[4]),
                             ston(shape_vec[5]));

        }

        obj.setPose(0.0, 0.0, 0.0); //default

        simulation[objects[i]] = obj;
    }

    //load POSES:

    //std::vector<std::string> poses = eclipse->query("getPoseList");
    std::vector<std::string> poses = swipl_interface::list2vector(swipl_interface::functor2vector(swi->query("findall(pose(P,X,Y,W),pose(P,X,Y,W),L)"))[3]);

    //for each pose
    for (int i = 0; i < poses.size(); i++) {
        //get the elements of the predicate
        std::vector<std::string> pose_vec = swipl_interface::functor2vector(poses[i]);
        //store the pose into the pose_map to be retreived during planning
        pose_map[pose_vec[1]] = Pose3d(ston(pose_vec[2]),
                                       ston(pose_vec[3]),
                                       ston(pose_vec[4]) );

        std::cout<<"POSE: "<<ston(pose_vec[2])<<" "<<ston(pose_vec[3])<<" "<<ston(pose_vec[4])<<std::endl;
    }

    //load map - if exists -
    //std::vector<std::string> map_vec = eclipse->query("getMap");
    std::vector<std::string> map_vec = swipl_interface::list2vector(swipl_interface::functor2vector(swi->query("map(M)"))[1]);
    if (map_vec.size() > 0) {
        map_file = map_vec[0];
        map_file.erase(std::remove(map_file.begin(), map_file.end(), '"'), map_file.end());
    } else {
        map_file = "";
        std::cout << "WARNING: map was not present in the PROLOG file" << std::endl;
    }


    //load the GOAL state:
    std::unordered_map < std::string, bool> tS_goal;

    //std::vector<std::string> goal_from_ltm = eclipse->query("getGoalState");
    std::vector<std::string> goal_from_ltm = swipl_interface::list2vector(swipl_interface::functor2vector(swi->query("goal_state(GS)"))[1]);

    if (goal_from_ltm.size() > 0) {
        for (int i = 0; i < goal_from_ltm.size(); i++)
            tS_goal[goal_from_ltm[i]] = true;

        Pose3d bS_goal(0.0, 0.0, 0.0); //this is more or less useless
        //Pose3d bS_goal(-2,0,0); //pose p2
        S_goal = State(tS_goal, bS_goal);
    } else
        std::cout << "WARNING: goal-state was not present in the PROLOG file" << std::endl;


    //load the INITIAL state:
    std::unordered_map < std::string, bool> tS_init;

    //std::vector<std::string> init_from_ltm = eclipse->query("getInitialState");
    std::vector<std::string> init_from_ltm = swipl_interface::list2vector(swipl_interface::functor2vector(swi->query("initial_state(GS)"))[1]);

    if (init_from_ltm.size() > 0) {
        for (int i = 0; i < init_from_ltm.size(); i++)
            tS_init[init_from_ltm[i]] = true;

        Pose3d bS_init(x, y, yaw);
        S_init = State(tS_init, bS_init);
        S_init.cout(false);
    } else
        std::cout << "WARNING: initial-state was not present in the PROLOG file" << std::endl;

}



//find the pose of an object, it can recursively search if an object is "on" another one!
//  NOTE: this function is DOMAIN-DEPENDENT !!!
Pose3d TM_RRTplanner::find_pose(std::unordered_map<std::string, bool>& var, std::string toFind) {
    std::vector<std::string> varVector;

    if (toFind == "")
        std::cout << "WHAT?! null target detected!" << std::endl;
    //if the target is a pose (e.g. p1) get it from the pose_map
    if (pose_map.find(toFind) != pose_map.end())
        return pose_map[toFind];

    //otherwise, it is an object, than find in the state the pose in which the object lies
    //for each predicate of the state
    for (auto it : var) {
        //parse the predicate 
        varVector = swipl_interface::functor2vector(it.first);
        //if the predicate is a on(X,Y) one (which means that and object X is "on" an object/pose Y)
        if (it.second && varVector[0] == "on" && varVector[1] == toFind) {
            if (pose_map.find(toFind) != pose_map.end()) //this if should be USELESS
                return pose_map[toFind];
            else
            //find the pose of the underlying object
            return find_pose(var, varVector[2]);
        }
    }
    //std::cout<<"TM_PLANNER WARNING: pose of "<<toFind<<" not found!"<<std::endl;
    return Pose3d(0, 0, 0);
}

bool TM_RRTplanner::is_pose_reached(Pose3d p1, Pose3d p2) {
    double lin_dist = sqrt(((p1.x - p2.x)*(p1.x - p2.x))+((p1.y - p2.y)*(p1.y - p2.y)));
    double ang_dist = std::abs(p1.w - p2.w);

    //if(lin_dist<=0.02 && ang_dist<=0.5)
    if (lin_dist <= 0.01 && ang_dist <= 0.5) //simulation
        return true;

    return false;
}

std::string TM_RRTplanner::get_object_loacation(State &state, std::string object_to_find) {
    std::vector<std::string> v;
    for (auto it : state.var) {
        if (it.second) {
            v = swipl_interface::functor2vector(it.first);
            if (v[0] == "on" && v[1] == object_to_find)
                return v[2];
        }
    }
    //std::cout<<ansi::red<<"LOCATION OF OBJECT '"<<object_to_find<<"' NOT FOUND!!"<<ansi::end<<std::endl;
    return "";
}

std::vector<Point3d> TM_RRTplanner::compute_dynamic_obstacles(State &state) {

    std::vector<Point3d> out_obs;

    for (auto it : simulation) {

        if (cout_debug_on) std::cout << "searching for " << it.first << " pose" << std::endl;

        //carrying
        if (!state.var["carry(" + it.first + ")"]) {

            if (cout_debug_on) std::cout << it.first << " not carried" << std::endl;

            out_obs = add_object_to_obstacles(it.first, state, out_obs);
        }
    }

    return out_obs;
}

std::vector<Point3d> TM_RRTplanner::compute_dynamic_obstacles(State &state, std::vector<Point3d> &obs) {

    std::vector<Point3d> out_obs = obs;

    for (auto it : simulation) {
        //carrying
        if (!state.var["carry(" + it.first + ")"])
            out_obs = add_object_to_obstacles(it.first, state, out_obs);
    }

    return out_obs;
}

//WARNING: this function is context specific!!!

std::vector<Point3d> TM_RRTplanner::add_object_to_obstacles(std::string obj_name, State &state) {

    std::string loc = get_object_loacation(state, obj_name);

    Pose3d pose;

    if (loc == "") //the cart has been taken
        pose = state.pose; //the cart is in the same position of the robot
    else
        pose = pose_map[loc]; //take the position of the cart

    if (cout_debug_on)
        std::cout << "pose of object " << obj_name << " is " << loc << ": (" << pose.x << "," << pose.y << "," << pose.w << ")" << std::endl;


    //simulation[obj_name].setPose(pose);
    //std::vector<Point3d> new_obstacles = simulation[obj_name].samplePoints(0.01);

    std::vector<Point3d> new_obstacles = simulation[obj_name].samplePointsFromPose(pose, 0.1);

    //update the obstacles
    return new_obstacles;

}

//WARNING: this function is context specific!!!

std::vector<Point3d> TM_RRTplanner::add_object_to_obstacles(std::string obj_name, State &state, std::vector<Point3d> &obs) {

    std::string loc = get_object_loacation(state, obj_name);

    Pose3d pose;

    if (loc == "") //the cart has been taken
        pose = state.pose; //the cart is in the same position of the robot
    else
        pose = pose_map[loc]; //take the position of the cart

    if (cout_debug_on)
        std::cout << "pose of object " << obj_name << " is " << loc << ": (" << pose.x << "," << pose.y << "," << pose.w << ")" << std::endl;

    //simulation[obj_name].setPose(pose);
    //std::vector<Point3d> new_obstacles = simulation[obj_name].samplePoints(0.01);

    std::vector<Point3d> new_obstacles = simulation[obj_name].samplePointsFromPose(pose, 0.1);

    new_obstacles.insert(new_obstacles.end(), obs.begin(), obs.end());

    //update the obstacles
    return new_obstacles;

}

std::vector<Point3d> TM_RRTplanner::remove_carrying_from_obstacles(State &state, std::vector<Point3d> &obs, double consis_margin) {

    std::vector<std::string> v;
    std::string carrying_obj = "";

    for (auto it : state.var) {
        if (it.second) {
            v = swipl_interface::functor2vector(it.first);
            if (v[0] == "carry") {
                carrying_obj = v[1];
                break;
            }
        }
    }

    if (carrying_obj == "")
        return obs;

    return remove_object_from_obstacles(carrying_obj, state, obs, consis_margin);
}

//WARNING: this function is context specific!!!

std::vector<Point3d> TM_RRTplanner::remove_object_from_obstacles(std::string obj_name, State &state, std::vector<Point3d> &obs, double consis_margin) {

    std::vector<Point3d> new_obstacles;

    //cart is called 'c' 
    std::string loc = get_object_loacation(state, obj_name);

    Pose3d pose;

    if (loc == "") //the cart has been taken
        pose = state.pose; //the cart is in the same position of the robot
    else
        pose = pose_map[loc]; //take the position of the cart

    //rotatate reference to be parallel with the pose
    double x_cart = (pose.x * std::cos(dtor(-pose.w)) - pose.y * std::sin(dtor(-pose.w)));
    double y_cart = (pose.y * std::cos(dtor(-pose.w)) + pose.x * std::sin(dtor(-pose.w)));

    double x_obs;
    double y_obs;

    double dx, dy;
    //for each obstacle
    for (auto i = 0; i < obs.size(); i++) {
        //rotate the reference of the robot to be parallel with the pose
        x_obs = (obs[i].x * std::cos(dtor(-pose.w)) - obs[i].y * std::sin(dtor(-pose.w)));
        y_obs = (obs[i].y * std::cos(dtor(-pose.w)) + obs[i].x * std::sin(dtor(-pose.w)));
        //simply compute the distance as the difference between the components and the robot semi-lenghts
        dx = std::abs(x_obs - x_cart) - ROBOT_CARRYING_LENGTH / 2;
        dy = std::abs(y_obs - y_cart) - ROBOT_CARRYING_WIDTH / 2;
        //if distances are smaller than the margin
        if ((dx > consis_margin) || (dy > consis_margin)) {
            //the pose is outside ..so add to the new obstacles
            new_obstacles.push_back(obs[i]);
        }
    }
    //update the obstacles
    return new_obstacles;

}

/** check if the pose is consistent adapting robot len considering the state
 *
 * @param state: the state to be checked pose and variables
 * @return 0 if the pose is colliding, otherwise the distance with the closer
 *      obstacle is returned.
 * 
 * //WARNING: this function is context specific!!!
 */
double TM_RRTplanner::is_pose_consistent(State &state, double consis_margin) {

    std::vector<Point3d> updated_obstacles = compute_dynamic_obstacles(state);

    updated_obstacles.insert(updated_obstacles.end(), obstacles.begin(), obstacles.end());

    //        std::cout<<"add_obs: "<<tester.toc()<<std::endl;

    if (state.var["carrying"])
        //return is_pose_consistent(state.pose, ROBOT_CARRYING_LENGTH, ROBOT_CARRYING_WIDTH, remove_carrying_from_obstacles(S_init, updated_obstacles, 0.2), consis_margin);
        return is_pose_consistent(state.pose, ROBOT_CARRYING_LENGTH, ROBOT_CARRYING_WIDTH, updated_obstacles, consis_margin);
    else
        return is_pose_consistent(state.pose, ROBOT_DEFAULT_LENGTH, ROBOT_DEFAULT_WIDTH, updated_obstacles, consis_margin);

}

/** check if the pose is consistent adapting robot len considering the state
 *
 * @param state: the state to be checked pose and variables
 * @param obstacles: the obstacles you can find in the current state
 * @return 0 if the pose is colliding, otherwise the distance with the closer
 *      obstacle is returned.
 * 
 * //WARNING: this function is context specific!!!
 */
double TM_RRTplanner::is_pose_consistent(State &state, std::vector<Point3d> &updated_obstacles, double consis_margin) {

    if (state.var["carrying"])
        //return is_pose_consistent(state.pose, ROBOT_CARRYING_LENGTH, ROBOT_CARRYING_WIDTH, remove_carrying_from_obstacles(S_init, updated_obstacles, 0.2), consis_margin);
        return is_pose_consistent(state.pose, ROBOT_CARRYING_LENGTH, ROBOT_CARRYING_WIDTH, updated_obstacles, consis_margin);
    else
        return is_pose_consistent(state.pose, ROBOT_DEFAULT_LENGTH, ROBOT_DEFAULT_WIDTH, updated_obstacles, consis_margin);

}

/** check if the pose is consistent (ie. not colliding obstacles)
 *
 * @param pose: the pose to be checked (x,y,yaw)
 * @param r_len, r_wid: robot length and robot width respectively
 * @return 0 if the pose is colliding, otherwise the distance with the closer
 *      obstacle is returned.
 */
double TM_RRTplanner::is_pose_consistent(Pose3d pose, double r_len, double r_wid, std::vector<Point3d> &obs, double consis_margin) {

    //rotatate reference to be parallel with the pose
    double x_rob = (pose.x * std::cos(dtor(-pose.w)) - pose.y * std::sin(dtor(-pose.w)));
    double y_rob = (pose.y * std::cos(dtor(-pose.w)) + pose.x * std::sin(dtor(-pose.w)));

    double x_obs;
    double y_obs;
    double min_obstacles_distance = 100;
    //margin of safety (added to the dimensions of the robot)

    double dx, dy;

    //std::cout<<"number of obstacles: "<<obs.size()<<std::endl;

    //for each obstacle
    for (auto i = 0; i < obs.size(); i++) {
        //rotate the reference of the robot to be parallel with the pose
        x_obs = (obs[i].x * std::cos(dtor(-pose.w)) - obs[i].y * std::sin(dtor(-pose.w)));
        y_obs = (obs[i].y * std::cos(dtor(-pose.w)) + obs[i].x * std::sin(dtor(-pose.w)));
        //simply compute the distance as the difference between the components and the robot semi-lenghts
        dx = std::abs(x_obs - x_rob) - r_len / 2;
        dy = std::abs(y_obs - y_rob) - r_wid / 2;
        //if distances are smaller than the margin
        if ((dx <= consis_margin) && (dy <= consis_margin)) {
            //the pose is colliding
            min_obstacles_distance = 0;
            break;
        }
        //find the closer obstacle
        if (min_obstacles_distance > std::sqrt((dx * dx) + (dy * dy))) {
            min_obstacles_distance = std::sqrt((dx * dx) + (dy * dy));
        }
    }

    //return the distance
    return min_obstacles_distance;
}

void TM_RRTplanner::draw_map() {
    draw_map(S_init.pose, S_goal.pose);
}

void TM_RRTplanner::draw_map(Pose3d pose, Pose3d goal) {
    //clear last map
    image = cv::Mat(RRT_DRAW_MAP_SIZE, RRT_DRAW_MAP_SIZE, CV_8UC3, cv::Scalar(0, 0, 0)); //1cm is 1px -> 1m is 100px

    double x_obs;
    double y_obs;

    draw_robot_pose(0, 0, pose.w, ROBOT_DEFAULT_LENGTH, ROBOT_DEFAULT_WIDTH, cv::Scalar(255, 150, 150)); //robot pose

    // draw_robot_pose(2,0,0,cv::Scalar(150,150,150)); //fix pose
    // draw_robot_pose(0,2,90,cv::Scalar(150,150,150)); //fix pose
    // draw_robot_pose(-2,0,180,cv::Scalar(150,150,150)); //fix pose
    // draw_robot_pose(0,-2,-90,cv::Scalar(150,150,150)); //fix pose

    //for each laser point (obstacle or not)
    for (auto i = 0; i < laser_points.size(); i++) {
        //transform it in pose-frame
        x_obs = laser_points[i].x - pose.x;
        y_obs = laser_points[i].y - pose.y;
        //plot the point as yellow circle
        cv::circle(image, cv::Point2f(x_obs * 100 + RRT_DRAW_MAP_SIZE / 2, -y_obs * 100 + RRT_DRAW_MAP_SIZE / 2), 4, cv::Scalar(0, 255, 255), 2, 8, 0);
    }

    //for each obstacle
    for (auto i = 0; i < obstacles.size(); i++) {
        //transform it in pose-frame
        x_obs = obstacles[i].x - pose.x;
        y_obs = obstacles[i].y - pose.y;
        //fill the yellow circle with red, it is an obstacle!
        cv::circle(image, cv::Point2f(x_obs * 100 + RRT_DRAW_MAP_SIZE / 2, -y_obs * 100 + RRT_DRAW_MAP_SIZE / 2), 2, cv::Scalar(0, 0, 255), 2, 8, 0);
    }

    //trasform the goal in pose-frame
    double x_g = goal.x - pose.x;
    double y_g = goal.y - pose.y;
    //plot the goal as a green cicrcle
    cv::circle(image, cv::Point2f(x_g * 100 + RRT_DRAW_MAP_SIZE / 2, -y_g * 100 + RRT_DRAW_MAP_SIZE / 2), 2, cv::Scalar(0, 255, 0), 2, 8, 0);

    draw_robot_pose(x_g, y_g, goal.w, ROBOT_DEFAULT_LENGTH, ROBOT_DEFAULT_WIDTH, cv::Scalar(0, 255, 0)); //goal pose

}

void TM_RRTplanner::draw_points(Pose3d pose, std::vector<Point3d> p, cv::Scalar col) {

    double x_p;
    double y_p;

    //for each point
    for (auto i = 0; i < p.size(); i++) {
        //fill the yellow circle with red, it is an obstacle!
        cv::circle(image, cv::Point2f((p[i].x - pose.x)*100 + RRT_DRAW_MAP_SIZE / 2, -(p[i].y - pose.y)*100 + RRT_DRAW_MAP_SIZE / 2), 2, col, 2, 8, 0);
    }
}

void TM_RRTplanner::draw_points(std::vector<Point3d> p, cv::Scalar col, double m2px, double offset) {

    double x_p;
    double y_p;

    //for each point
    for (auto i = 0; i < p.size(); i++) {
        //fill the yellow circle with red, it is an obstacle!
        cv::circle(image, cv::Point2f(p[i].x * m2px + offset, -p[i].y * m2px + offset), 2, col, 2, 8, 0);
    }
}


void TM_RRTplanner::draw_plan(Pose3d pose, int start_index) {

    double x_p;
    double y_p;
    //for each step of the plan
    for (auto i = start_index; i < plan.size(); i++) {
        //std::cout<<"\t"<<i<<": "<<plan[i].pos.x<<", "<<plan[i].pos.y<<", "<<plan[i].pos.w<<std::endl;
        //plot the plan in pose-frame
        x_p = plan[i].state.pose.x - pose.x;
        y_p = plan[i].state.pose.y - pose.y;

        draw_robot_pose(x_p, y_p, plan[i].state.pose.w, plan[i].state, cv::Scalar(255, 0, 0)); //step pose

    }
    // //show the image with CV
    // cv::imshow("motionPlanner",image);
    // cv::waitKey(3);
}

void TM_RRTplanner::draw_plan(std::vector<PlanStepTM> in_plan, Pose3d pose) {

    double x_p;
    double y_p;
    //for each step of the plan
    for (auto i = 1; i < in_plan.size(); i++) {
        //std::cout<<"\t"<<i<<": "<<plan[i].pos.x<<", "<<plan[i].pos.y<<", "<<plan[i].pos.w<<std::endl;
        //plot the plan in pose-frame
        x_p = in_plan[i].state.pose.x - pose.x;
        y_p = in_plan[i].state.pose.y - pose.y;

        draw_robot_pose(x_p, y_p, in_plan[i].state.pose.w, in_plan[i].state, cv::Scalar(255, 0, 0) ); //step pose

    }
    // //show the image with CV
    // cv::imshow("motionPlanner",image);
    // cv::waitKey(3);
}

void TM_RRTplanner::draw_robot_pose(double px, double py, double pw, double r_len, double r_wid, cv::Scalar color, double m2px, double offset) {
    //get the rect representing the robot pose
    cv::RotatedRect rRect = cv::RotatedRect(cv::Point2f(px * m2px + offset, -py * m2px + offset), cv::Size2f(r_len*m2px, r_wid * m2px), -pw);
    cv::Point2f vertices[4];
    //plot the rect
    rRect.points(vertices);
    for (int i = 0; i < 4; i++)
        cv::line(image, vertices[i], vertices[(i + 1) % 4], color, 2, 8, 0);

    //cv::arrowedLine(image, cv::Point2f(px*m2px +offset, -py*m2px +offset), cv::Point2f(px*m2px +offset +std::cos(-pw)*10, -py*m2px +offset +std::sin(-pw)*10), color, 2, 8, 0, 0.2);
    //plot a circle corresponding to the robot frontal part
    cv::circle(image, cv::Point2f(vertices[2].x + (vertices[3].x - vertices[2].x) / 2, vertices[2].y + (vertices[3].y - vertices[2].y) / 2), 3, color, 2, 8, 0);
}

void TM_RRTplanner::draw_robot_pose(Pose3d p, double r_len, double r_wid, cv::Scalar color, double m2px, double offset) {
    draw_robot_pose(p.x, p.y, p.w, r_len, r_wid, color, m2px, offset);
}

void TM_RRTplanner::draw_robot_pose(State state, cv::Scalar color, double m2px, double offset) {
    if (state.var["carrying"])
        draw_robot_pose(state.pose.x, state.pose.y, state.pose.w, ROBOT_CARRYING_LENGTH, ROBOT_CARRYING_WIDTH, color, m2px, offset);
    else
        draw_robot_pose(state.pose.x, state.pose.y, state.pose.w, ROBOT_DEFAULT_LENGTH, ROBOT_DEFAULT_WIDTH, color, m2px, offset);
}

void TM_RRTplanner::draw_robot_pose(double px, double py, double pw, State state, cv::Scalar color, double m2px, double offset) {
    if (state.var["carrying"])
        draw_robot_pose(px, py, pw, ROBOT_CARRYING_LENGTH, ROBOT_CARRYING_WIDTH, color, m2px, offset);
    else
        draw_robot_pose(px, py, pw, ROBOT_DEFAULT_LENGTH, ROBOT_DEFAULT_WIDTH, color, m2px, offset);
}

double TM_RRTplanner::symmetric_difference(std::unordered_map<std::string, bool>& v1, std::unordered_map<std::string, bool>& v2) {
    //v1 is always the small one

    double d = 0;

    std::unordered_map < std::string, bool>::const_iterator fnd;
    //copy the small one (time optimization)
    std::unordered_map < std::string, bool> app = v1;

    //to speed-up the algo, always copy the small one
    for (auto it : v2) {
        fnd = app.find(it.first);
        if (fnd != app.end()) {
            if (fnd->second != it.second) {
                d++;
            }
            app.erase(fnd);
        } else if (it.second != false) {
            d++;
        }
    }
    //now check the rest of the app elements
    for (auto it : app) {
        if (v2[ it.first ] != it.second)
            d++;
    }

    return d;
}

void TM_RRTplanner::update_obstacles_from_LIDAR() {
    //get obstacles in world-frame
    //obstacles = robot2world_frame(WMV.get<std::vector < Point3d >> ("laser.obstacles"));
    //get laser readings in world-frame
    //laser_points = robot2world_frame(WMV.get<std::vector < Point3d >> ("laser.points"));
}

void TM_RRTplanner::update_obstacles_from_FILE(std::string path_to_map) {
    cv::Mat map_img;

    if (path_to_map == ""){
        map_img = cv::imread(path_to_directory + map_file, cv::IMREAD_GRAYSCALE);
        std::cout<<"map-file: "<<path_to_directory << map_file<<std::endl;
    }
    else{
        map_img = cv::imread(path_to_directory + path_to_map, cv::IMREAD_GRAYSCALE);
    }

    //double px2meters = 0.01; //this means 1 px -> 0.01 meters
    double px2meters = 0.1; //this means 1 px -> 0.1 meters
    //double px2meters = 0.05; //this means 1 px -> 0.05 meters

    double iy, jx;

    obstacles.clear();

    double y_off = map_img.cols / 2;
    double x_off = map_img.rows / 2;

    for (int j = 0; j < map_img.rows; j++) {
        jx = j;
        for (int i = 0; i < map_img.cols; i++) {
            iy = i;
            if (map_img.at<uchar>(j, i) < 128) { //black!
                //remember that the origin of the map is in the middle of the image!
                obstacles.push_back(Point3d(-(jx - x_off) * px2meters, -(iy - y_off) * px2meters, 0)); //now it should be in map frame!
            }
        }
    }

}



//transform a vector of points from robot to world frame

std::vector<Point3d> TM_RRTplanner::robot2world_frame(std::vector<Point3d> vec) {
    double app_x;
    double app_y;
    //for each point
    for (auto i = 0; i < vec.size(); i++) {
        //transform it in world frame (x,y,yaw are robot coordinates)
        app_x = x + (vec[i].x * std::cos(dtor(yaw)) - vec[i].y * std::sin(dtor(yaw)));
        app_y = y + (vec[i].y * std::cos(dtor(yaw)) + vec[i].x * std::sin(dtor(yaw)));

        vec[i].x = app_x;
        vec[i].y = app_y;
    }
    //return the vector
    return vec;
}


//see planner.cpp file for its decalration
//std::vector<PlanStepTM> TM_RRTplanner::plan_rrt_simple(State& start, State& goal, std::vector<PlanStepTM>& rrt_plan, double rrt_timeout, double horizon, double rrt_step) {


//OPTIMIZATION of RRT: find shortcuts bypassing noisy paths

std::vector<PlanStepTM> TM_RRTplanner::random_shortcut_plan(std::vector<PlanStepTM>& in_plan, double rs_timeout) {
    double stop_distance = 0.01;

    std::vector<PlanStepTM> out_plan;

    //if plan is too small
    if (in_plan.size() < 2)
        //do nothing
        return in_plan;

    out_plan = in_plan;
    //start the clock
    double elapsed_secs = 0;
    clock_t end_clock;
    clock_t begin_clock = clock();
    //while u have time
    while (elapsed_secs < rs_timeout) {

        //probability to optimize the initial part of the plan
        double P_optimize_init = 0.3;
        //flip a float-coint between 0 and 1
        double coin = fRand(0, 1);

        //get two random steps of the plan
        int r_start, r_stop;

        //if initial part is selected
        if (coin <= P_optimize_init) {
            r_start = 0;
            r_stop = iRand(r_start + 1, out_plan.size() - 1);
        } else {
            r_start = iRand(0, out_plan.size() - 2);
            r_stop = iRand(r_start + 1, out_plan.size() - 1);
        }

        // std::cout<<"[shortcut-id] start: "<<r_start<<" stop: "<<r_stop<<std::endl;

        State s_curr = out_plan[r_start].state;
        State s_stop = out_plan[r_stop].state;
        State s_goal = out_plan[out_plan.size() - 1].state;

        std::vector<PlanStepTM> shortcut;
        shortcut.clear();

        bool shortcut_broken = false;
        bool shortcut_found = false;
        //while the shortcut is not found and not broken (shortcut failed)
        while (!shortcut_broken && !shortcut_found) {
            //make a step from the last element of the shortcut to the stop pose
            StepMT step_new = action_in_direction(s_curr, s_stop, s_goal);
            State s_new;
            s_new.pose = estimate_new_pose(s_curr.pose, step_new.motion);
            s_new.var = s_curr.var;
            //if the shortcut is too big, we assume it is going nowhere
            if (shortcut.size() + 1 > 100) {
                //the current shortcut is broken
                shortcut_broken = true;
                continue;
            }
            //if the new step is free from obstacles
            double obs_dist = is_pose_consistent(s_new);
            if (obs_dist > 0) {
                //and does not reduce the distance from obstacles
                if (obs_dist < out_plan[out_plan.size() - 1].min_obst) {
                    shortcut_broken = true;
                    continue;
                }
                //check if the stop pose is reached
                if (distance(s_new.pose, s_stop.pose) <= 0.01 && distance(s_new.var, s_stop.var) == 0) {
                    //check if the shortcut is too small
                    if (shortcut.size() <= 1) {
                        //in case, it is useless
                        shortcut_broken = true;
                        continue;
                    }
                    //otherwise we've found a shortcut!
                    shortcut_found = true;
                    std::cout << ansi::yellow << "shortcut of size " << shortcut.size() << "found!" << ansi::end << std::endl;
                    //insert the shortcut into the plan, cutting off the old points
                    out_plan.erase(out_plan.begin() + r_start + 1, out_plan.begin() + r_stop);
                    for (auto i = 0; i < shortcut.size(); i++) {
                        out_plan.insert(out_plan.begin() + r_start + i + 1, shortcut[i]);
                    }

                }                    //oth. add the pose to the shortcut
                else {

                    // std::cout<<"p_new:  "<<p_new.x<<","<<p_new.y<<","<<p_new.w<<std::endl;
                    // std::cout<<"p_stop: "<<p_stop.x<<","<<p_stop.y<<","<<p_stop.w<<std::endl;
                    // std::cout<<"pos: "<<ps.pos.x<<","<<ps.pos.y<<","<<ps.pos.w<<std::endl;
                    // std::cout<<"act: "<<ps.act.fs<<","<<ps.act.ls<<","<<ps.act.ts<<std::endl;
                    // std::cout<<"shortcut-push"<<std::endl;
                    //compute the closest obstacle from the shortcut
                    if (shortcut.size() > 0 && shortcut[shortcut.size() - 1].min_obst < obs_dist) {
                        obs_dist = shortcut[shortcut.size() - 1].min_obst;
                    }
                    //insert the step
                    PlanStepTM ps(s_new, step_new, 0, obs_dist);
                    shortcut.push_back(ps);

                    //move forward
                    s_curr = s_new;
                }
            }                //oth. we collide an obstacle
            else {
                //the shortcut is brocken
                shortcut_broken = true;
            }
        }
        //update the clock
        end_clock = clock();
        elapsed_secs = double(end_clock - begin_clock) / CLOCKS_PER_SEC;

    }

    //return the new (hopefully better) plan
    return out_plan;
}

//get the plan-step that is closer to the current robot pose

int TM_RRTplanner::get_nearest_plan_step() {
    int near_index = -1;
    double near_lin_distance = 10;
    double near_ang_distance = 100;
    //get the robot pose
    Pose3d rob_pos(x, y, yaw);

    if (plan.size() == 0)
        return -1;

    //State rob_state(plan[curr_plan_step].state.var, rob_pos);
    State rob_state(S_init); //this should be called current ...not init!
    double near_state_distance = 10000;

    //for each plan step
    for (auto i = curr_plan_step; i < plan.size(); i++) {
        if (distance(plan[i].state, rob_state) < near_state_distance) {
            near_state_distance = distance(plan[i].state, rob_state);
            near_index = i;
        }

        /*
                //get the distance with the robot (considering the rotation)
                double lin_dist = distance2d(rob_pos,plan[i].state.pose);
                double ang_dist = std::abs(rob_pos.w - plan[i].state.pose.w);
                //if i is closer
                if(lin_dist<near_lin_distance){
                    //save i
                    near_lin_distance = lin_dist;
                    near_ang_distance = ang_dist;
                    near_index = i;
                }
                //oth. if i is the "same" position with better angle
                else if(std::abs(lin_dist - near_lin_distance)<0.01 && ang_dist<near_ang_distance){
                    //save i
                    near_lin_distance = lin_dist;
                    near_ang_distance = ang_dist;
                    near_index = i;
                }
         */
    }
    return near_index;
}

bool TM_RRTplanner::check_plan() {
    bool consis = true;
    //if already have a plan
    for (auto i = curr_plan_step; i < plan.size(); i++) {
        //if plan is no more consistent
        //if(!is_pose_consistent(plan[i].state.pose,robot_length,robot_width,0.0)){
        if (!is_pose_consistent(plan[i].state, 0.0)) {
            //skip the rest of the plan
            std::cout << "plan collide at step " << i << std::endl;
            consis = false;
            break;
        }
    }
    return consis;
}

std::unordered_map<std::string, bool> TM_RRTplanner::apply_to_state(Task& t, State& s) {
    return apply_to_state(t, s.var);
}

std::unordered_map<std::string, bool> TM_RRTplanner::apply_to_state(Task& t, std::unordered_map<std::string, bool>& v) {
    std::vector<std::string> cond = t.post_conditions;
    std::unordered_map < std::string, bool> new_v = v;

    for (auto i = 0; i < cond.size(); i++) {
        if (cond[i].substr(0, 1) == "-")
            new_v[ cond[i].substr(1, cond[i].size()) ] = false;
        else
            new_v[ cond[i] ] = true;
    }
    return new_v;

}

bool TM_RRTplanner::is_applicable(Task& t, State& s) {
    return is_applicable(t, s.var);
}

bool TM_RRTplanner::is_applicable(Task& t, std::unordered_map<std::string, bool>& inv) {
    std::vector<std::string> cond = t.pre_conditions;
    std::unordered_map < std::string, bool> v = inv;

    bool consis = true;
    int i = 0;

    while (i < cond.size() && consis) {
        if ((cond[i].substr(0, 1) == "-" && v[ cond[i].substr(1, cond[i].size()) ] == true) ||
                cond[i].substr(0, 1) != "-" && v[ cond[i] ] == false) {
            consis = false;
        }
        i++;
    }
    return consis;
}

StepMT TM_RRTplanner::action_in_direction(State& s1, State& s2, State& g) {
    //top-side
    Task best_task = task_in_direction(s1.var, s2.var);
    //bottom-side
    Step3d best_action = step_in_direction(s1.pose, s2.pose, g.pose);
    return StepMT(best_action, best_task);

}


/**** **** **** TASK IN DIRECTION SELECTOR **** **** ****/

Task TM_RRTplanner::task_in_direction(State& s1, State& s2, TaskSelector mode){
    switch( mode ){
        case TaskSelector::BEST:
            return task_in_direction_BEST(s1.var, s2.var);
        case TaskSelector::UNIFORM:
            return task_in_direction_UNIF(s1.var, s2.var);
        case TaskSelector::MONTECARLO:
            return task_in_direction_MC(s1.var, s2.var);
    }
}

Task TM_RRTplanner::task_in_direction(std::unordered_map<std::string,bool>& v1, std::unordered_map<std::string,bool>& v2, TaskSelector mode){
    switch( mode ){
        case TaskSelector::BEST:
            return task_in_direction_BEST(v1, v2);
        case TaskSelector::UNIFORM:
            return task_in_direction_UNIF(v1, v2);
        case TaskSelector::MONTECARLO:
            return task_in_direction_MC(v1, v2);
    }
}


/**** **** **** TASK IN DIRECTION BEST **** **** ****/
    
Task TM_RRTplanner::task_in_direction_BEST(State& s1, State& s2){
    return task_in_direction_BEST(s1.var, s2.var);
}

Task TM_RRTplanner::task_in_direction_BEST(std::unordered_map<std::string,bool>& v1, std::unordered_map<std::string,bool>& v2){

    //NOTE: also in this case we can find different tasks with the same distance
    std::vector<Task> vec_best_task;

    Task best_task;
    std::unordered_map<std::string,bool> best_var_state;
    double best_distance = -1;
    for(auto i=0; i<task_set.size(); i++){
        //check if task_set[i].precond are satisfied by v1
        if( is_applicable(task_set[i], v1) ){
            //compute the new_var_state applying task_set[i].postcond to s1.var
            std::unordered_map<std::string,bool> new_var_state = apply_to_state(task_set[i], v1);
            //compute the distance betweend new_var_state and s2.var with the best distance
            double new_distance = distance(new_var_state, v2);
            //save the task with the best distance
            if( best_distance<0 || new_distance < best_distance ){
                best_task = task_set[i];
                best_var_state = new_var_state;
                best_distance = new_distance;

                vec_best_task.clear();
                vec_best_task.push_back(task_set[i]);
            }
            else if(new_distance == best_distance){
                vec_best_task.push_back(task_set[i]);
            }
        }
    }

    //select randomly the best task among those with the same (lowest) distance
    int last_index = vec_best_task.size()-1;

    if(last_index>0)
        best_task = vec_best_task[ iRand(0,last_index) ];
    else
        best_task = vec_best_task[0];

//        //PLOT
//        if(going_to_goal && best_distance == 0){
//            last_task_selected = true;
//            std::cout<<ansi::yellow<<"last task selected: "<<best_task.name<<ansi::end<<std::endl;
//        }
//        else
//            last_task_selected = false;

    return best_task;

}

/**** **** **** TASK IN DIRECTION UNIFORM **** **** ****/

Task TM_RRTplanner::task_in_direction_UNIF(State& s1, State& s2){
    return task_in_direction_UNIF(s1.var, s2.var);
}

Task TM_RRTplanner::task_in_direction_UNIF(std::unordered_map<std::string,bool>& v1, std::unordered_map<std::string,bool>& v2){

    //select randomly the best task among those with the same (lowest) distance
    int last_index = task_set.size()-1;

    Task best_task = task_set[ iRand(0,last_index) ];

    int limit = 5;

    while(!is_applicable(best_task, v1) && limit>0){
        best_task = task_set[ iRand(0,last_index) ];
        limit--;
    }

    if(limit<=0)
        best_task.name = ""; //return no task!

    return best_task;
}

/**** **** **** TASK IN DIRECTION MONTECARLO **** **** ****/

Task TM_RRTplanner::task_in_direction_MC(State& s1, State& s2){
    return task_in_direction_MC(s1.var, s2.var);
}

Task TM_RRTplanner::task_in_direction_MC(std::unordered_map<std::string,bool>& v1, std::unordered_map<std::string,bool>& v2){

    //select randomly the best task among those with the same (lowest) distance
    int last_index = task_set.size()-1;

    Task best_task;
    best_task.name = ""; //return no task!

    Task candidate_task;

    int limit = task_set.size()/3;


    std::unordered_map<std::string,bool> best_var_state;
    double best_distance = -1;

    while(limit>0){
        //sample a random task
        candidate_task = task_set[ iRand(0,last_index) ];
        //check if task_set[i].precond are satisfied by v1
        if( is_applicable(candidate_task, v1) ){
            //compute the new_var_state applying task_set[i].postcond to s1.var
            std::unordered_map<std::string,bool> new_var_state = apply_to_state(candidate_task, v1);
            //compute the distance betweend new_var_state and s2.var with the best distance
            double new_distance = distance(new_var_state, v2);
            //save the task with the best distance
            if( new_distance == 0 ){
                best_task = candidate_task;
                break;
            }
            else if( best_distance<0 || new_distance < best_distance || (new_distance == best_distance && iRand(0,1)==1) ){
                best_task = candidate_task;
                best_var_state = new_var_state;
                best_distance = new_distance;

            }
        }
        limit--;
    }

    return best_task;
}

Step3d TM_RRTplanner::step_in_direction(State& s1, State& s2, State& g) {
    return step_in_direction(s1.pose, s2.pose, g.pose);
}

//steer functions (RRT)

Step3d TM_RRTplanner::step_in_direction(Pose3d p1, Pose3d p2, Pose3d goal) {

    //step sizes
    double max_fs_step = 0.1; //m
    double max_ls_step = 0.05; //m
    double max_ts_step = 5; //deg

    //compute goal distances
    double goal_dist = sqrt(((goal.x - p1.x)*(goal.x - p1.x))+((goal.y - p1.y)*(goal.y - p1.y)));
    double goal_ang_dist = std::abs(goal.w - p1.w);
    // if( goal_dist < 0.5 ){
    //     //decrease the step size!
    //     max_fs_step = max_fs_step * (goal_dist/0.5);
    //     max_ls_step = max_ls_step * (goal_dist/0.5);
    //     max_ts_step = max_ts_step * (goal_dist/0.5);
    // }


    //compute new step
    deg180 delta_yaw;
    deg180 delta_yaw_rob;

    Step3d s; //new step

    //if the goal is close
    if (goal_dist <= 0.05) {
        //go directly to the goal
        delta_yaw = rtod(std::atan2((goal.y - p1.y), (goal.x - p1.x)));
        delta_yaw_rob = delta_yaw - p1.w;

        //s.fs = goal_dist * std::cos(dtor(delta_yaw)); //std::cos(dtor(delta_yaw_rob)); //changed 02/12/2021 to correct the trajectory overshooting 
        //s.ls = goal_dist * std::sin(dtor(delta_yaw)); //std::sin(dtor(delta_yaw_rob));

        s.fs = goal_dist * std::cos(dtor(delta_yaw_rob)); //changed back 15/12/2021
        s.ls = goal_dist * std::sin(dtor(delta_yaw_rob));

        s.ts = goal.w - p1.w;
        s.ts = std::abs(s.ts) > max_ts_step ? sgn(s.ts) * max_ts_step : (double) s.ts;
    }        //oth. move in p2 direction
    else {

        delta_yaw = rtod(std::atan2((p2.y - p1.y), (p2.x - p1.x)));
        delta_yaw_rob = delta_yaw - p1.w;

        double dist_lin = sqrt(((p2.x - p1.x)*(p2.x - p1.x))+((p2.y - p1.y)*(p2.y - p1.y)));

        s.fs = dist_lin * std::cos(dtor(delta_yaw_rob)); //std::cos(dtor(delta_yaw)); //std::cos(dtor(delta_yaw_rob));
        s.fs = std::abs(s.fs) > max_fs_step ? sgn(s.fs) * max_fs_step : s.fs;

        s.ls = dist_lin * std::sin(dtor(delta_yaw_rob)); //std::sin(dtor(delta_yaw)); //std::sin(dtor(delta_yaw_rob));
        s.ls = std::abs(s.ls) > max_ls_step ? sgn(s.ls) * max_ls_step : s.ls;

        s.ts = p2.w - p1.w;
        s.ts = std::abs(s.ts) > max_ts_step ? sgn(s.ts) * max_ts_step : (double) s.ts;
    }

    return s;
}


//states have to be the same vars!!

std::vector< PlanStepTM > TM_RRTplanner::path_in_direction(State& s_start, State& s_stop, Pose3d& p_goal, double stop_distance) {

    //std::vector<PlanStepTM> path;
    partial_path.clear();

    State s_curr = s_start;

    double min_obs = -1;

    //compute the dynamic obstacles for the whole path (optimization)
    //  (they only depends on vars, which are constant)
    std::vector<Point3d> state_obstacles = compute_dynamic_obstacles(s_start);
    state_obstacles.insert(state_obstacles.end(), obstacles.begin(), obstacles.end());

    //while the path is not too long (ie. is not looping)
    while (partial_path.size() < 100 * stop_distance) {
        //make a step from the last element of the path to the stop pose
        PlanStepTM step_new;
        step_new.act.motion = step_in_direction(s_curr.pose, s_stop.pose, p_goal);
        //State s_new;
        step_new.state.pose = estimate_new_pose(s_curr.pose, step_new.act.motion);
        step_new.state.var = s_curr.var;
        //if the new step is free from obstacles
        //double obs_dist = is_pose_consistent(step_new.state); //this was computing dynamic obstacles every time
        double obs_dist = is_pose_consistent(step_new.state, state_obstacles);
        if (obs_dist > 0) {

            if (min_obs < 0 || obs_dist < min_obs) {
                min_obs = obs_dist;
                step_new.min_obst = min_obs;
            }

            //check if the stop pose is reached
            if (is_pose_reached(step_new.state.pose, s_stop.pose) ||
                    //or the goal pose is reached
                    is_pose_reached(step_new.state.pose, p_goal) ||
                    //or the path limit is reached
                    distance(step_new.state.pose, s_start.pose) >= stop_distance) {
                //we've found a path!
                partial_path.push_back(step_new);

                return partial_path; //complete path
            }                
            //oth. add the pose to the path
            else {
                partial_path.push_back(step_new);
                //move forward
                s_curr = step_new.state;
            }
        }            
        //oth. we collide an obstacle
        else {
            //return path; //partial path
            return std::vector<PlanStepTM>(); //empty path
        }
    }

    path_looping_count++;

    return std::vector<PlanStepTM>(); //empty path
}

Pose3d TM_RRTplanner::estimate_new_pose(Pose3d p, Step3d s) {
    Pose3d p_new; //compute new position
    double freq = 1; //sec

    double delta_fs_y = std::sin(dtor(p.w)) * s.fs*freq;
    double delta_fs_x = std::cos(dtor(p.w)) * s.fs*freq;
    double delta_ls_y = std::sin(dtor(p.w + 90)) * s.ls*freq;
    double delta_ls_x = std::cos(dtor(p.w + 90)) * s.ls*freq;

    p_new.x = p.x + delta_fs_x + delta_ls_x;
    p_new.y = p.y + delta_fs_y + delta_ls_y;

    p_new.w = p.w + s.ts*freq;

    return p_new;
}

void TM_RRTplanner::rviz_plot_plan(std::vector< PlanStepTM > vec) {

    geometry_msgs::Point p;

    visualization_msgs::Marker markers;

    markers.header.frame_id = "/base_link"; //"/map";
    markers.header.stamp = ros::Time::now();
    markers.ns = "plan";
    markers.action = visualization_msgs::Marker::ADD;
    markers.pose.orientation.w = 1.0;

    markers.id = 0;

    markers.type = visualization_msgs::Marker::POINTS;

    // POINTS markers use x and y scale for width/height respectively
    markers.scale.x = 0.1;
    markers.scale.y = 0.1;

    // Points are
    markers.color.g = 1.0f;
    markers.color.a = 1.0;

    for (auto i = 0; i < vec.size(); i++) {

        p.x = vec[i].state.pose.x;
        p.y = vec[i].state.pose.y;
        p.z = 0.0;


        markers.points.push_back(p);

    }

    marker_pub.publish(markers);
}

std::vector< PlanStepTM > TM_RRTplanner::transform_plan(std::vector< PlanStepTM > vec, std::string in_frame, std::string out_frame) {

    std::vector< PlanStepTM > res;

    geometry_msgs::PointStamped p_in, p_out;


    //laser_listener.waitForTransform(in_frame, out_frame, ros::now, ros::Duration(0.5));


    //wait for transform with map (NOTE: the mutex must be unlocked!)
    if (!laser_listener.waitForTransform(in_frame, out_frame, ros::Time(0), ros::Duration(0.01))) {
        return res;
    }

    for (auto i = 0; i < vec.size(); i++) {

        p_in.header.frame_id = in_frame;
        p_in.point.x = vec[i].state.pose.x;
        p_in.point.y = vec[i].state.pose.y;
        p_in.point.z = vec[i].state.pose.w;

        laser_listener.transformPoint(out_frame, ros::Time(0), p_in, in_frame, p_out);

        State app(vec[i].state.var, Pose3d(p_out.point.x, p_out.point.y, p_out.point.z));

        res.push_back(PlanStepTM(app, vec[i].act, vec[i].cost, vec[i].min_obst));
    }
    return res;
}

void TM_RRTplanner::rviz_image_plan() {

    //std::cout<<"publishing image_plan"<<std::endl;

    sensor_msgs::ImagePtr msg = cv_bridge::CvImage(std_msgs::Header(), "bgr8", image).toImageMsg();

    image_pub.publish(msg);
}


void TM_RRTplanner::ros_publish_plan() {

    std_msgs::String msg;

    std::stringstream ss;
    ss<<"( ";
    for(auto i=1; i<plan.size(); i++){
        ss << plan[i].toString();
        if(i != plan.size()-1)
            ss<<",";
    }
    ss<<")";
    
    msg.data = ss.str();
    plan_pub.publish(msg);

    std::cout<<"PLAN:"<<std::endl;
    std::cout<<ss.str()<<std::endl;
}


void TM_RRTplanner::draw_state(State s, cv::Scalar color, double m2px, double offset){
    
    double r_len = ROBOT_DEFAULT_LENGTH;
    double r_wid = ROBOT_DEFAULT_WIDTH;
    
    if (s.var["carrying"]){
        r_len = ROBOT_CARRYING_LENGTH;
        r_wid = ROBOT_CARRYING_WIDTH;
    }
    
    //get the rect representing the robot pose
    cv::RotatedRect rRect = cv::RotatedRect(cv::Point2f(s.pose.x * m2px + offset, -s.pose.y * m2px + offset), cv::Size2f(r_len*m2px, r_wid * m2px), -s.pose.w);
    cv::Point2f vertices[4];
    //plot the rect
    rRect.points(vertices);
    for (int i = 0; i < 4; i++)
        cv::line(image, vertices[i], vertices[(i + 1) % 4], color, 2, 8, 0);
    
    //plot a circle corresponding to the robot frontal part
    cv::circle(image, cv::Point2f(vertices[2].x + (vertices[3].x - vertices[2].x) / 2, vertices[2].y + (vertices[3].y - vertices[2].y) / 2), 3, color, 2, 8, 0);
    
    cv::putText(image, s.toString(true,true) , vertices[0], cv::FONT_HERSHEY_DUPLEX, 0.65, color, 2);
}

void TM_RRTplanner::draw_path(std::vector<PlanStepTM> path, cv::Scalar color, double m2px, double offset){
    
    double r_len;
    double r_wid;

    
    for(auto p=0; p<path.size(); p++){
        r_len = ROBOT_DEFAULT_LENGTH;
        r_wid = ROBOT_DEFAULT_WIDTH;

        if (path[p].state.var["carrying"]){
            r_len = ROBOT_CARRYING_LENGTH;
            r_wid = ROBOT_CARRYING_WIDTH;
        }
    
        //get the rect representing the robot pose
        cv::RotatedRect rRect = cv::RotatedRect(cv::Point2f(path[p].state.pose.x * m2px + offset, -path[p].state.pose.y * m2px + offset), cv::Size2f(r_len*m2px, r_wid * m2px), -path[p].state.pose.w);
        cv::Point2f vertices[4];
        //plot the rect
        rRect.points(vertices);
        for (int i = 0; i < 4; i++)
            cv::line(image, vertices[i], vertices[(i + 1) % 4], color, 2, 8, 0);
        
        //plot a circle corresponding to the robot frontal part
        cv::circle(image, cv::Point2f(vertices[2].x + (vertices[3].x - vertices[2].x) / 2, vertices[2].y + (vertices[3].y - vertices[2].y) / 2), 3, color, 2, 8, 0);
    }
}

void TM_RRTplanner::draw_sample(State s_near, State s_new, State s_sample){
    
    draw_state(s_near, cv::Scalar(255, 100, 0));
    
    draw_state(s_new, cv::Scalar(255, 200, 0));
    
    draw_state(s_sample, cv::Scalar(200, 200, 200));
    
    rviz_image_plan();
}

void TM_RRTplanner::draw_sample(State s_near, State s_sample){
    
    draw_state(s_near, cv::Scalar(255, 100, 0));
    
    draw_state(s_sample, cv::Scalar(200, 200, 200));
    
    rviz_image_plan();
}

void TM_RRTplanner::draw_tree(RRTree rrt, cv::Scalar color, double m2px, double offset){
    State cs, ps;
    
    draw_map();
    
    for(auto i=1; i<rrt.size(); i++){
        cs = rrt.node[i];
        ps = rrt.node[ rrt.parent[i] ];
        
        cv::circle(image, cv::Point2f(cs.pose.x * m2px + offset, -cs.pose.y * m2px + offset), 3, color, 2, 8, 0);
        
        cv::line(image, cv::Point2f(cs.pose.x * m2px + offset, -cs.pose.y * m2px + offset), cv::Point2f(ps.pose.x * m2px + offset, -ps.pose.y * m2px + offset), color, 2, 8, 0);
    }
}


// on(c2,p3) draw_cluster_poses(RRTree t, RRTcluster c, "on(c2,p3)")
//returns the id of the cluster

int TM_RRTplanner::draw_cluster_poses(RRTree t, RRTcluster c, std::vector<std::string> var) {
    int id = -1;
    for (int i = 0; i < c.nodes.size(); i++) {

        int fnd = 0;

        for (auto k = 0; k < var.size(); k++) {
            //for each element of the cluster var_map
            for (auto it : c.vars[i]) {
                //if the element is present in the var_map
                if (it.first == var[k] && it.second) {
                    //the cluster is different
                    fnd++;
                    id = i;
                    break;
                }
            }
        }

        if (fnd != var.size())
            continue;

        //std::cout << "DRAWING CLUSTER -------------> " << id << std::endl;
        std::stringstream ss;
        ss << "cluster: "<<id;
        cv::putText(image, ss.str() ,cv::Point(15, 30) , cv::FONT_HERSHEY_DUPLEX, 0.65, cv::Scalar(0, 143, 143), 2);
        for (int j = 0; j < c.nodes[id].size(); j++) {
            draw_robot_pose(t.node[c.nodes[id][j]], cv::Scalar(255, 0, 200));
        }

        draw_points(S_init.pose, compute_dynamic_obstacles(t.node[c.nodes[id][0]]));

        return id;
    }
}

void TM_RRTplanner::draw_cluster_poses(RRTree t, RRTcluster c, int id) {
    //std::cout << "DRAWING CLUSTER -------------> " << id << std::endl;
    std::stringstream ss;
    ss << "cluster: "<<id;
    cv::putText(image, ss.str() ,cv::Point(15, 30) , cv::FONT_HERSHEY_DUPLEX, 0.65, cv::Scalar(0, 143, 143), 2);
    for (int i = 0; i < c.nodes[id].size(); i++) {
        draw_robot_pose(t.node[c.nodes[id][i]], cv::Scalar(255, 0, 200));
    }
}

void TM_RRTplanner::cout_cluster(RRTcluster c, State _goal, bool truthset_only) {

    //quick and dirty map in alphabetic order
    std::map < std::string, bool, std::greater < std::string>> ordered;

    for (int i = 0; i < c.vars.size(); i++) {
        //ordering
        ordered.clear();
        for (auto it : c.vars[i]) {
            ordered[it.first] = it.second;
        }

        std::cout << "CLUSTER" << i << "( ( ";
        for (auto it = ordered.begin(); it != ordered.end(); it++) {
            if (!truthset_only)
                std::cout << it->first << ":" << it->second << " ";
            else if (it->second)
                std::cout << it->first << " ";
        }
        std::cout << ")" << std::endl;

        std::cout<<"\t number of samples : "<<c.sampled[i]<<std::endl;
        std::cout<<"\t rejected due task: "<<c.rejected_due_task[i]<<std::endl;
        std::cout<<"\t rejected due path: "<<c.rejected_due_path[i]<<std::endl;
        std::cout<<"\t nodes in the cluster: "<<c.nodes[i].size()<<std::endl;
        std::cout<<"\t distance from goal: "<<distance(c.vars[i],_goal.var)<<std::endl;
	
    }
}

//lpots only id, number of nodes and distance

void TM_RRTplanner::cout_less_cluster(RRTcluster c, State _goal) {
    for (int i = 0; i < c.vars.size(); i++) {
        std::cout << "CLUSTER" << i << ": " << c.nodes[i].size() << std::endl;
        std::cout << "\t   distance from goal: " << distance(c.vars[i], _goal.var) << std::endl;
    }
}

double TM_RRTplanner::median(std::vector<double> scores) {
    double median;
    size_t size = scores.size();

    std::sort(scores.begin(), scores.end());

    if (size % 2 == 0) {
        median = (scores[size / 2 - 1] + scores[size / 2]) / 2;
    } else {
        median = scores[size / 2];
    }

    return median;
}

