using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace rBergomi
{
    //public interface IreturnPrice
    //{
    //    double price();
    //}

    public class VIXOption
    {
        public enum OptionType { Call = 1, Put = 0 }
        public double maturity { get; set; }
        public double strike { get; set; }
        public OptionType type{ get; set; }

        public VIXOption(OptionType type_,double maturity_,double strike_)
        {
            type = type_;
            maturity = maturity_;
            strike = strike_;
        }

    }
}
