using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace MathNet.Numerics.SystemModels
{
    /// <summary>
    /// This interface indicates, that a system model can be simulated over one dimension (t) given n inputs (u) and m states (x)
    /// by solving dx(u,x,t) = f(u,x,t) and y(u,x,t) = f(u,x,t)
    /// </summary>
    public interface ISimulateable
    {
        /// <summary>
        /// method to call once before simulation starts
        /// </summary>
        /// <param name="startTime">the starting time for this simulation</param>
        /// <param name="p"></param>
        void Initialize(double startTime, params object[] p);

        /// <summary>
        /// Model name, can be null or empty
        /// </summary>
        string ModelName { get; set; }

        /// <summary>
        /// Input signal names must be same size as
        /// </summary>
        string[] InputNames { get; set; }
        string[] OutputNames { get; set; }
        string[] StateNames { get; set; }

        /// <summary>
        /// the number of states in the model
        /// </summary>
        int nStates { get; }

        /// <summary>
        /// the number of inputs for the model
        /// </summary>
        int nInputs { get; }

        /// <summary>
        /// the number of outputs from the model
        /// </summary>
        int nOutputs { get; }

        /// <summary>
        /// Called first when simulating, should calculate y based on x, u and t
        /// </summary>
        /// <param name="u">input vector for this model at time t</param>
        /// <param name="x">state vector for this model at time t</param>
        /// <param name="t">time now t</param>
        /// <returns>the outputs of this model as calculated by f(u,x,t)</returns>
        Vector<double> CalcOutputs(Vector<double> u, Vector<double> x, double t);

        /// <summary>
        /// Called second, when simulating, should calculate x_k+1 (state for the next step) based on x, u and t
        /// </summary>
        /// <param name="u">input vector for this model at time t</param>
        /// <param name="x">state vector for this model at time t</param>
        /// <param name="t">time now t</param>
        /// <returns>the derivative of this models dynamics as calculated by f(u,x,t)</returns>
        Vector<double> CalcDiscreteStates(Vector<double> u, Vector<double> x, double t);

        /// <summary>
        /// Called third, when simulating, should calculate dx based on x, u and t
        /// </summary>
        /// <param name="u">input vector for this model at time t</param>
        /// <param name="x">state vector for this model at time t</param>
        /// <param name="t">time now t</param>
        /// <returns>the derivative of this models dynamics as calculated by f(u,x,t)</returns>
        Vector<double> CalcDerivatives(Vector<double> u, Vector<double> x, double t);

        void Finalize(params object[] p);

    }
}
