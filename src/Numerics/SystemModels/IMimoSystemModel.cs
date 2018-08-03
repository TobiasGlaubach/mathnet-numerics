using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;

namespace MathNet.Numerics.SystemModels
{
    public interface IMimoSystemModel
    {
        string[] InputNames { get; }
        string[] OutputNames { get; }
        bool IsStable();
        IEnumerable<double[]> CalcResponse(IEnumerable<double[]> uVecs);
        IEnumerable<Vector<double>> CalcResponse(IEnumerable<Vector<double>> uVecs);

        IEnumerable<Vector<double>> Impulse(int nSteps);
        IEnumerable<Vector<double>> Impulse();
        IEnumerable<Complex[]> Bode(int nPoints = 100);
        IEnumerable<Complex[]> Bode(int nPoints, out double[] omega_vec);
        IEnumerable<Complex[]> Bode(double[] omega_vec);
        IEnumerable<Complex[]> GetPoles();
        IEnumerable<Complex[]> GetZeros();
    }
}
