#include <iostream>
#include <eigen3/Eigen/Dense>


//convert radiant to degree
inline double rtod(double r){
    return r * 180.0/ M_PI;
}

//convert degree to radiant
inline double dtor(double r){
    return r * M_PI / 180.0;
}


class deg180{
public:
    deg180(){
        this->ang = 0;
    }
    deg180(double a){
        this->ang = ( a>180 ? a-360 : (a<=-180 ? a+360 : a) );
    }

    // Overload (double) operator to return the angle value.
    inline operator double() { return this->ang; }

    // Overload + operator to add two ang180.
    inline double operator+(double a)
    {
        //std::cout<<"deg180: "<<this->ang<<" + "<<a<<" = "<<this->ang<<std::endl;
        double res = this->ang + a;
        res = ( res>180 ? res-360 : (res<=-180 ? res+360 : res) );

        return res;
    }

    // Overload - operator to subtradt two ang180.
    inline double operator-(double a)
    {
        //std::cout<<"deg180: "<<this->ang<<" - "<<a<<" = "<<this->ang<<std::endl;
        double res = this->ang - a;
        res = ( res>180 ? res-360 : (res<=-180 ? res+360 : res) );

        return res;
    }

    // Overload = operator to set from double.
    inline double operator=(double a)
    {
        this->ang = ( a>180 ? a-360 : (a<=-180 ? a+360 : a) );
        return this->ang;
    }
    inline double rad(){
        return ( this->ang * M_PI / 180.0 );
    }
private:
    double ang;
};

class Point3d{
public:
    Point3d(){
        x = 0;
        y = 0;
        z = 0;
    };
    Point3d(double nx, double ny, double nz){
        x = nx;
        y = ny;
        z = nz;
    };
    Point3d(double nx, double ny){
        x = nx;
        y = ny;
        z = 0; //std::nan;
    };
    double x; //m
    double y; //m
    double z; //m
};

class Pose3d{
public:
    Pose3d(){
        x = 0;
        y = 0;
        w = 0;
    };
    Pose3d(double nx, double ny, deg180 nw){
        x = nx;
        y = ny;
        w = nw;
    };
    Pose3d(double nx, double ny){
        x = nx;
        y = ny;
        w = 0; //std::nan;
    };
    // Overload = operator to set from Pose3d.
    inline Pose3d operator=(Pose3d a)
    {
        this->x = a.x;
        this->y = a.y;
        this->w = a.w;
        return a;
    }
    double x; //m
    double y; //m
    //double w; //deg
    deg180 w; //deg
};


template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

    

///SIMULATION PART THIS NEEDS EIGEN!!

/**
 * Turn the Pose3d into homogeneous transform matrix to facilitate transforms
 * 
 * @param p: Pose3d to be transformed in homogeneous-matrix
 * @return the 3x3 homogeneous transform matrix in Eigen
 */
inline Eigen::MatrixXd pose2HT(Pose3d p){ //don't know why but the coordinates are inverted!
    Eigen::MatrixXd t(4,4);
    Eigen::MatrixXd r(4,4);
    double a = -p.w.rad(); //dtor(p.w);
    // Matrix starts from (0,0)
    r << std::cos(a), -std::sin(a), 0.0, 0.0,
         std::sin(a), std::cos(a),  0.0, 0.0,
            0.0,          0.0,      1.0, 0.0,
            0.0,          0.0,      0.0, 1.0;
    
    t <<    1.0,          0.0,      0.0, -p.x, 
            0.0,          1.0,      0.0, -p.y,
            0.0,          0.0,      1.0, 0.0,
            0.0,          0.0,      0.0, 1.0;
    
    return r*t;
}

#define SHAPE_MIN_SAMPLING_STEP 0.001

class Shape2d{
public:
    Shape2d(){
        pose.x = 0.0;
        pose.y = 0.0;
        pose.w = 0.0;
        T_os = pose2HT(pose);
        T_so = T_os.inverse();
    }
    void setPose(Pose3d p){
        pose = p;
        T_os = pose2HT(pose);
        T_so = T_os.inverse();
    }
    void setPose(double nx, double ny, double nw){
        pose.x = nx;
        pose.y = ny;
        pose.w = nw;
        T_os = pose2HT(pose);
        T_so = T_os.inverse();
    }
    inline Pose3d getPose(){
        return pose;
    }
    inline Eigen::MatrixXd getTso(){
        return T_so;
    }
    inline Eigen::MatrixXd getTos(){
        return T_os;
    }
    inline std::vector<Point3d> getPoints(double stp){
        return samples;
    }
    virtual std::vector<Point3d> samplePoints(double stp) = 0;
//private:
    Eigen::MatrixXd T_os; //transform from outside to shape_pose
    Eigen::MatrixXd T_so; //transform from shape_pose to outside
    Pose3d pose;
    
    double sampling_step; //in meters!
    std::vector<Point3d> samples;
};


//the base shape2d is a square
class ShapeRect : public Shape2d{
public:
    ShapeRect(){
        pose.x = 0.0;
        pose.y = 0.0;
        pose.w = 0.0;
        T_os = pose2HT(pose);
        T_so = T_os.inverse();
    }
    ShapeRect(Pose3d p, double xl, double yl, double stp = 0.1){
        pose = p;
        T_os = pose2HT(pose);
        T_so = T_os.inverse();
        x_len = xl;
        y_len = yl;
        
        samplePoints(stp);
    }
    ShapeRect(double x, double y, deg180 w, double xl, double yl, double stp = 0.1){
        pose.x = x;
        pose.y = y;
        pose.w = w;
        T_os = pose2HT(pose);
        T_so = T_os.inverse();
        x_len = xl;
        y_len = yl;
        
        samplePoints(stp);
    }
    
    std::vector<Point3d> samplePoints(double stp){
        
        stp = stp < SHAPE_MIN_SAMPLING_STEP ? SHAPE_MIN_SAMPLING_STEP : stp;
        
//        if(stp == sampling_step)
//            return samples;
//        sampling_step = stp;
        
        int x_stp = x_len/stp < 1 ? 1 : x_len/stp;
        int y_stp = y_len/stp < 1 ? 1 : y_len/stp;
        
        samples.clear();
        
        Eigen::VectorXd v(4);
        
        //sample x-lines
        for(int i=0; i<x_stp; i++){
            //transform from shape to out frame
            v(0) = -x_len/2.0 + stp * ( (double) i );
            v(1) = y_len/2.0;
            v(2) = 0.0;
            v(3) = 1.0;
            v = T_so * v;
            //add to the samples
            samples.push_back( Point3d( v(0), v(1) ) );
        
            if(y_len <= 0.0)
                continue;
            
            //transform from shape to out frame
            v(0) = x_len/2.0 - stp * ( (double) i );
            v(1) = -y_len/2.0;
            v(2) = 0.0;
            v(3) = 1.0;
            v = T_so * v;
            //add to the samples
            samples.push_back( Point3d( v(0), v(1) ) );
        }
        
        //sample y-lines
        for(int i=0; i<y_stp; i++){
            //transform from shape to out frame
            v(0) = x_len/2.0;
            v(1) = y_len/2.0 - stp * ( (double) i );
            v(2) = 0.0;
            v(3) = 1.0;
            v = T_so * v;
            //add to the samples
            samples.push_back( Point3d( v(0), v(1) ) );
        
            if(x_len <= 0.0)
                continue;
            
            //transform from shape to out frame
            v(0) = -x_len/2.0;
            v(1) = -y_len/2.0 + stp * ( (double) i );
            v(2) = 0.0;
            v(3) = 1.0;
            v = T_so * v;
            //add to the samples
            samples.push_back( Point3d( v(0), v(1) ) );
        }
        
        return samples;
    }
    
    double x_len;
    double y_len;
};

class Object2d : public Shape2d{
public:
    Object2d(){
        pose.x = 0.0;
        pose.y = 0.0;
        pose.w = 0.0;
        T_os = pose2HT(pose);
        T_so = T_os.inverse();
    }
    Object2d(Pose3d p, double stp = 0.1){
        pose = p;
        T_os = pose2HT(pose);
        T_so = T_os.inverse();
        
        samplePoints(stp);
    }
    Object2d(double x, double y, deg180 w, double stp = 0.1){
        pose.x = x;
        pose.y = y;
        pose.w = w;
        T_os = pose2HT(pose);
        T_so = T_os.inverse();
        
        samplePoints(stp);
    }
    
    std::vector<Point3d> samplePoints(double stp){
        
        stp = stp < SHAPE_MIN_SAMPLING_STEP ? SHAPE_MIN_SAMPLING_STEP : stp;
        
//        if(stp == sampling_step)
//            return samples;
//        sampling_step = stp;
        
        samples.clear();
        
        Eigen::VectorXd v(4);
        
        //sample x-lines
        for(int i=0; i<shapes.size(); i++){
            
            std::vector<Point3d> shape_samples = shapes[i].samplePoints(stp);
            
            for(int j=0; j<shape_samples.size(); j++){
                //transform from shape to out frame
                v(0) = shape_samples[j].x;
                v(1) = shape_samples[j].y;
                v(2) = 0.0;
                v(3) = 1.0;
                v = T_so * v;
                //add to the samples
                samples.push_back( Point3d( v(0), v(1) ) );
            }
        }
        
        return samples;
        
    }
    
    std::vector<Point3d> samplePointsFromPose(Pose3d p, double stp){
        
        stp = stp < SHAPE_MIN_SAMPLING_STEP ? SHAPE_MIN_SAMPLING_STEP : stp;
        
        Eigen::MatrixXd t_so(4,4);
        t_so = pose2HT(p).inverse();
        
        Eigen::VectorXd v(4);
        
        std::vector<Point3d> out_samples;
        
        //sample x-lines
        for(int i=0; i<shapes.size(); i++){
            
            std::vector<Point3d> shape_samples = shapes[i].samplePoints(stp);
            
            for(int j=0; j<shape_samples.size(); j++){
                //transform from shape to out frame
                v(0) = shape_samples[j].x;
                v(1) = shape_samples[j].y;
                v(2) = 0.0;
                v(3) = 1.0;
                v = t_so * v;
                //add to the samples
                out_samples.push_back( Point3d( v(0), v(1) ) );
            }
        }
        
        return out_samples;
        
    }
    
    inline void addShapeRect(Pose3d p, double xl, double yl, double stp = 0.1){
        shapes.push_back( ShapeRect(p,xl,yl,stp) );
    }
    inline void addShapeRect(double x, double y, deg180 w, double xl, double yl, double stp = 0.1){
        shapes.push_back( ShapeRect(x,y,w,xl,yl,stp) );
    }
    inline std::vector<ShapeRect> getShapes(){
        return shapes;
    }
private:
    std::vector<ShapeRect> shapes;
};

///END OF SIMULATION PART