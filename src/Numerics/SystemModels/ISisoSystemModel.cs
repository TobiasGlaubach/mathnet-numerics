using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;

namespace MathNet.Numerics.SystemModels
{
    public interface ISisoSystemModel
    {
        string Name { get; }
        bool IsStable();
        IEnumerable<double> CalcResponse(IEnumerable<double> x);
        double[] CalcResponse(double[] x);
        double[] Impulse(int nSteps);
        double[] Impulse();
        Complex[] Bode(int nPoints = 100);
        Complex[] Bode(int nPoints, out double[] omega_vec);
        Complex[] Bode(double[] omega_vec);
        Complex[] GetPoles();
        Complex[] GetZeros();

    }
}
