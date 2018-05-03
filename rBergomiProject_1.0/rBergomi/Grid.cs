using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace rBergomi
{
    public class Grid
    {
        double timeStep_;// delta _t
        double t_min_; //t_0
        int timeNmbrStep_; // N nombre d'echantillage temporel

        public Grid(double t_min, double t_max, int timeNmbrStep)
        {
            t_min_ = t_min;
            timeNmbrStep_ = timeNmbrStep;
            timeStep_ = (double) (t_max - t_min) / timeNmbrStep;
        }

        //functions
        public  double t(int i)
        {
            if (i<0||i > timeNmbrStep_) throw new IndexOutOfRangeException();
            return t_min_ + i * timeStep_;
        } 

        //getters
        public int get_timeNmbrStep() { return timeNmbrStep_; }
        public double get_Step() { return timeStep_; }
    }
}
