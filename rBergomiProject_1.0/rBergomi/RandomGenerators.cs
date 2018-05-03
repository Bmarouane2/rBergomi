using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Meta.Numerics.Statistics.Distributions;

namespace rBergomi
{
    public interface IUniformSimulator
    {
        double next();
    }
    
    public class UniformSimulator:IUniformSimulator
    {
        private readonly Random random_;

        public UniformSimulator(Random random)
        {
            random_ = random;
        }

        public double next()
        {
            return random_.NextDouble();
        }
    }

    public class GaussianSimulator
    {
        private readonly IUniformSimulator random_ = new UniformSimulator(new Random());
        private readonly NormalDistribution normalLaw_ = new NormalDistribution();

        public double Next()
        {
            return normalLaw_.InverseRightProbability(random_.next());
        }
    }
}
