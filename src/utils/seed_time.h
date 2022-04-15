#include <iostream>
#include <chrono>

/*
 *  seed::time is a wrapper for chrono, the most precise time-counter
 *      as far as I know.
 */

namespace seed {
namespace time {

//a point in the time
typedef std::chrono::high_resolution_clock::time_point t_point;

//elapsed seconds between two time points
inline double elapsed(t_point t0, t_point t1){
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    return time_span.count();
}

//get the current time point (like now() function)
inline t_point now(){
    return std::chrono::high_resolution_clock::now();
}

//stopwatch matlab-like implementing tic and toc
class Clock {
public:
    Clock(){
        t0 = seed::time::now();
    }

    t_point tic(){
        t0 = seed::time::now();
        return t0;
    }

    //returns elapsed time in seconds
    double toc(){
        t_point t1 = seed::time::now();
        return seed::time::elapsed(t0, t1);
    }

private:
    t_point t0;
};

//associate the time to a value of a generic type (template)
//  the time is set iff the value is set!
template<class T>
class TimedValue {
public:
    TimedValue(){
        time = seed::time::now();
        value = T();
    }
    TimedValue(T val){
        time = seed::time::now();
        value = val;
    }
    TimedValue(T val, t_point t){
        time = t;
        value = val;
    }
    void set(T val){
        time = seed::time::now();
        value = val;
    }
    T get(){
        return value;
    }
    t_point getTime(){
        return time;
    }
protected:
    t_point time;
    T value;
};

}
}
