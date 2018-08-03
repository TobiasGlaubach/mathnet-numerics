using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using MathNet.Numerics.SystemModels;

namespace MathNet.Numerics.SystemModels.LtiSystems
{
    /// <summary>
    /// a class for LTI state space models in the form
    /// x_k1 = A * x_k  + B * u
    /// y_k1 = C * x_k1 + D * u
    /// </summary>
    public class StateSpaceDiscrete : IMimoSystemModel, IDiscreteSystemModel
    {
        
        private Matrix<double> _A;

        /// <summary>
        /// State transition matrix
        /// </summary>
        public Matrix<double> A
        {
            get { return _A; }
            set { checkDimensions(value, B, C, D, value.ColumnCount, B.ColumnCount, C.ColumnCount); _A = value; }
        }

        private Matrix<double> _B;

        /// <summary>
        /// Input to state matrix
        /// </summary>
        public Matrix<double> B
        {
            get { return _B; }
            set { checkDimensions(A, value, C, D, A.ColumnCount, value.ColumnCount, C.ColumnCount); _B = value; }
        }

        private Matrix<double> _C;

        /// <summary>
        /// State to output matrix
        /// </summary>
        public Matrix<double> C
        {
            get { return _C; }
            set { checkDimensions(A, B, value, D, A.ColumnCount, B.ColumnCount, value.ColumnCount); _C = value; }
        }

        
        private Matrix<double> _D;

        /// <summary>
        /// Feedthrough Matrix (input to output matrix)
        /// </summary
        public Matrix<double> D
        {
            get { return _D; }
            set { checkDimensions(A, B, C, value, A.ColumnCount, B.ColumnCount, C.ColumnCount); _D = value; }
        }


        private Vector<double> _x_k1;

        /// <summary>
        /// current states
        /// </summary>
        public Vector<double> x_k1
        {
            get { return _x_k1; }
            set
            {
                if (value.Count != nStates)
                {
                    throw new ArgumentException("dimension given does not match the expected number of states (exp {nStates}, got: {value.Count}");
                }
                _x_k1 = value;
            }
        }

        /// <summary>
        /// Sampling time for this state space model
        /// </summary>
        public double Ts { get; set; }

        private string[] _StateNames;

        /// <summary>
        /// the names for this models states
        /// </summary>
        public string[] StateNames
        {
            get { return _StateNames; }
            set
            {
                if (value.Length != nStates)
                {
                    throw new ArgumentException("dimension given does not match the expected number of states (exp: {nStates}, got: {value.Count}");
                }

                _StateNames = value;
            }
        }

        private string[] _InputNames;

        /// <summary>
        /// the names for this models inputs
        /// </summary>
        public string[] InputNames
        {
            get { return _InputNames; }
            set
            {
                if (value.Length != nInputs)
                {
                    throw new ArgumentException("dimension given does not match the expected number of inputs (exp: {nInputs}, got: {value.Count}");
                }

                _InputNames = value;
            }
        }

        private string[] _OutputNames;

        /// <summary>
        /// the names for this models outputs
        /// </summary>
        public string[] OutputNames
        {
            get { return _OutputNames; }
            set
            {
                if (value.Length != nOutputs)
                {
                    throw new ArgumentException("dimension given does not match the expected number of outputs (exp: {nOutputs}, got: {value.Count}");
                }
                _OutputNames = value;
            }
        }

        /// <summary>
        /// gets the number of states from A.ColumnCount
        /// </summary>
        public int nStates { get { return A.ColumnCount; } }

        /// <summary>
        /// gets the number of inputs from B.ColumnCount
        /// </summary>
        public int nInputs { get { return B.ColumnCount; } }

        /// <summary>
        /// gets the number of outputs from C.ColumnCount
        /// </summary>
        public int nOutputs { get { return C.ColumnCount; } }


        public StateSpaceDiscrete()
        {

        }

        /// <summary>
        /// Constructs a new state space model based on another model (in a new memory location)
        /// </summary>
        /// <param name="G"></param>
        /// <param name="copyState">if set to true, x_k1 from G will be copied as well</param>
        public StateSpaceDiscrete(StateSpaceDiscrete G, bool copyState = false)
        {
            init(G.A, G.B, G.C, G.D, copyData: true);
            StateNames = G.StateNames?.Select(n => (string)n.Clone()).ToArray();
            InputNames = G.InputNames?.Select(n => (string)n.Clone()).ToArray();
            OutputNames = G.OutputNames?.Select(n => (string)n.Clone()).ToArray();
            Ts = G.Ts;
            if (copyState)
            {
                x_k1 = DenseVector.OfVector(G.x_k1);
            }
        }

        public StateSpaceDiscrete(int nStates, int nInputs = 1, int nOutputs = 1)
        {
            if (nStates <= 0)
                throw new ArgumentOutOfRangeException("nStates must be greater than zero (got: " + nStates + ")");

            if (nInputs <= 0)
                throw new ArgumentOutOfRangeException("nInputs must be greater than zero (got: " + nInputs + ")");

            if (nOutputs <= 0)
                throw new ArgumentOutOfRangeException("nOutputs must be greater than zero (got: " + nOutputs + ")");

            var a = DenseMatrix.CreateDiagonal(nStates, nStates, 1.0);
            var b = DenseMatrix.CreateDiagonal(nStates, nInputs, 1.0);
            var c = DenseMatrix.Create(nOutputs, nStates, 1.0);
            var d = new DenseMatrix(nOutputs, nInputs);
            init(a, b, c, d, copyData: false);
        }

        /// <summary>
        /// creates a new siso state space model with just one in- and output as well as N states
        /// </summary>
        /// <param name="A">state transition matrix</param>
        /// <param name="copyData">if true, a new block of memory will be allocated for each matrix given, if false, the matrices will be used directly (faster but may result in unforseen effects)</param>
        public StateSpaceDiscrete(Matrix<double> A, bool copyData=true)
        {
            if (A == null)
                throw new ArgumentNullException("A");

            var arr = new double[A.ColumnCount];
            for (int i = 0; i < arr.Length; i++)
                arr[i] = 1.0d;

            var B = DenseMatrix.Create(1, 1, 1.0);
            var C = DenseMatrix.Create(1, A.ColumnCount, 1.0);
            var D = DenseMatrix.Create(1, 1, 0.0);
            init(A, B, C, D, copyData);
        }

        /// <summary>
        /// creates a new state space model from given system matrices
        /// </summary>
        /// <param name="A">state transition matrix</param>
        /// <param name="B">input to state matrix</param>
        /// <param name="C">state to output matrix</param>
        /// <param name="D">direct feedthrough matrix (input to output matrix)</param>
        /// <param name="copyData">if true, a new block of memory will be allocated for each matrix given, if false, the matrices will be used directly (faster but may result in unforseen effects)</param>
        public StateSpaceDiscrete(Matrix<double> A, Matrix<double> B, Matrix<double> C, Matrix<double> D, bool copyData=true)
        {
            init(A, B, C, D, copyData);
        }

        private void init(Matrix<double> A, Matrix<double> B, Matrix<double> C, Matrix<double> D, bool copyData)
        {
            if (A == null)
                throw new ArgumentNullException("A");
            if (B == null)
                throw new ArgumentNullException("B");
            if (C == null)
                throw new ArgumentNullException("C");
            if (D == null)
                throw new ArgumentNullException("D");

            checkDimensions(A, B, C, D, A.ColumnCount, B.ColumnCount, C.ColumnCount);
            if (copyData)
            {
                _A = DenseMatrix.OfMatrix(A);
                _B = DenseMatrix.OfMatrix(B);
                _C = DenseMatrix.OfMatrix(C);
                _D = DenseMatrix.OfMatrix(D);
            }
            else
            {
                _A = A;
                _B = B;
                _C = C;
                _D = D;
            }

            x_k1 = new DenseVector(nStates);
        }


        private static void checkDimensions(Matrix<double> A, Matrix<double> B, Matrix<double> C, Matrix<double> D, int nStates, int nInputs, int nOutputs)
        {
            string str = "{0} Matrix {1} count must be equal to {2} (expected: {3}, got:{4})";

            if (A.ColumnCount != nStates)
                throw new ArgumentException(String.Format(str, "A", "column", "nStates", nStates, A.ColumnCount));

            if (A.RowCount != nStates)
                throw new ArgumentException(String.Format(str, "A", "row", "nStates", nStates, A.RowCount));

            if (B.ColumnCount != nInputs)
                throw new ArgumentException(String.Format(str, "B", "column", "nInputs", nInputs, B.ColumnCount));

            if (B.RowCount != nStates)
                throw new ArgumentException(String.Format(str, "B", "column", "nInputs", nStates, B.RowCount));

            if (A.ColumnCount != nStates)
                throw new ArgumentException(String.Format(str, "C", "column", "nStates", nStates, C.ColumnCount));

            if (A.RowCount != nOutputs)
                throw new ArgumentException(String.Format(str, "C", "row", "nOutputs", nOutputs, C.RowCount));

            if (D.ColumnCount != nInputs)
                throw new ArgumentException(String.Format(str, "D", "column", "nInputs", nInputs, D.ColumnCount));

            if (D.RowCount != nOutputs)
                throw new ArgumentException(String.Format(str, "D", "column", "nOutputs", nOutputs, D.RowCount));

        }

        #region Operators


        private static Matrix<double> combineDiag(Matrix<double> M1, Matrix<double> M2)
        {
            var M = new DenseMatrix(M1.RowCount + M2.RowCount, M1.ColumnCount + M2.ColumnCount);
            M.SetSubMatrix(0, 0, M1);
            M.SetSubMatrix(M1.RowCount, M1.ColumnCount, M2);

            return M;           
        }

        private static string[] makeStringNames(int size, string name = null)
        {
            var lst = new string[size];
            for (int i = 0; i < size; i++)
                lst[i] = name + i.ToString();
            return lst;
        }


        public void DropOutput(string outputToDrop)
        {
            if (OutputNames == null)
            {
                throw new NotSupportedException("You can't drop outputs by using an output name if your outputs are not named");
            }

            if (!OutputNames.Contains(outputToDrop))
            {
                throw new KeyNotFoundException("the output name you tried to drop ({outputToDrop}) was not found in the output names");
            }

            var C_tmp = C;
            var D_tmp = D;
            var outputNames = OutputNames.Where(n => !n.Equals(outputToDrop));

            int cnt = 0;
            foreach (var name in OutputNames)
            {
                if (name.Equals(outputToDrop))
                {
                    C_tmp = C.RemoveRow(cnt);
                    D_tmp = D.RemoveRow(cnt);
                }
                else
                    cnt++;
            }

            OutputNames = null;
            init(A, B, C_tmp, D_tmp, copyData: false);
            OutputNames = outputNames.ToArray();
        }

        public void DropOutput(int index)
        {
            if (index <= 0 || index >= nOutputs)
            {
                throw new ArgumentOutOfRangeException("index");
            }

            var C_tmp = C;
            var D_tmp = D;
            var outputNames = OutputNames.ToList();
            outputNames.RemoveAt(index);

            C_tmp = C.RemoveRow(index);
            D_tmp = D.RemoveRow(index);

            OutputNames = null;
            init(A, B, C_tmp, D_tmp, copyData: false);
            OutputNames = outputNames.ToArray();
        }

        /// <summary>
        /// Link (also called cascade, or series) two state space models to each other by mapping the output of the 
        /// first state space to the input of the second state space.
        /// Input and Output names must match if present. Internal states will be disregarded
        /// </summary>
        /// <param name="G1">the first model whichs outputs to link to the second model</param>
        /// <param name="G2">the second model whichs outputs to link to the second model</param>
        /// <param name="shrinkNames"></param>
        /// <see cref="https://math.stackexchange.com/questions/2201067/cascade-of-state-space-models-for-linear-systems"/>
        /// <returns>a new independent state space model</returns>
        public static StateSpaceDiscrete Link(StateSpaceDiscrete G1, StateSpaceDiscrete G2, bool shrinkNames = true)
        {
            if ((G1.OutputNames == null) != (G2.InputNames == null))
                throw new NotSupportedException("Linking a non named outputs with a named inputs or vice versa is not supported, either drop all names first, or make both models equal in regard with names");


            if (G1.OutputNames != null && G2.InputNames != null)
            {
                if(!G1.OutputNames.SequenceEqual(G2.InputNames))
                    throw new NotSupportedException("In order to link two named state space models with each other, the first models output names must equal the second models input names. Drop all names in both models if you want to ignore this.");
            }
            else
            {
                if (G1.nOutputs != G2.nInputs)
                    throw new NotSupportedException("In order to link two state space models with each other, the first models number of output must match the second models inputs");
            }

            if (G1.Ts != G2.Ts)
            {
                throw new ArgumentException("Link not possible, since the two models you tried to link have different sampling times ({G1.Ts) vs {G2.Ts}");
            }

            // T.Glaubach: this could be done more efficiently, but this implementation improves readability and makes it quite easy to grasp whats going on
            var A1 = G1.A;
            var A2 = G2.A;
            var B1 = G1.B;
            var B2 = G2.B;
            var C1 = G1.C;
            var C2 = G2.C;
            var D1 = G1.D;
            var D2 = G2.D;

            var A21 = B2 * C1;
            var A22 = A2;
            var A11 = A1;
            var A12 = new DenseMatrix(A11.RowCount, A22.ColumnCount);


        }

        /// <summary>
        /// Combine (mux two parallel models) two state space models into one by diagonally joining their in- and output Matrix as well as state matrix in one matrix if necessary.
        /// In- Output and Statenames not taken into account. 
        /// </summary>
        /// <param name="G1">first state space model</param>
        /// <param name="G2">second spate space model</param>
        /// <returns></returns>
        public static StateSpaceDiscrete Combine(StateSpaceDiscrete G1, StateSpaceDiscrete G2)
        {

            List<string> sNames = null;
            List<string> iNames = null;
            List<string> oNames = null;


            if ((G1.StateNames == null) != (G2.StateNames == null))
                throw new NotSupportedException("Combining a non named state with a named state is not supported");
            if ((G1.InputNames == null) != (G2.InputNames == null))
                throw new NotSupportedException("Combining a non named state with a named state is not supported");
            if ((G1.OutputNames == null) != (G2.OutputNames == null))
                throw new NotSupportedException("Combining a non named state with a named state is not supported");

            var forceNamed = false;
            if ((G1.OutputNames != null) || (G1.StateNames != null) || (G1.InputNames != null) ||
                (G2.OutputNames != null) || (G2.StateNames != null) || (G2.InputNames != null))
            {
                forceNamed = true;
            }

            if (forceNamed)
            {
                sNames = new List<string>();
                iNames = new List<string>();
                oNames = new List<string>();

                sNames.AddRange(G1.StateNames ?? makeStringNames(G1.nStates, "G1_"));
                sNames.AddRange(G2.StateNames ?? makeStringNames(G2.nStates, "G2_"));

                iNames.AddRange(G1.InputNames ?? makeStringNames(G1.nInputs, "G1_"));
                iNames.AddRange(G2.InputNames ?? makeStringNames(G2.nInputs, "G2_"));

                oNames.AddRange(G1.OutputNames ?? makeStringNames(G1.nOutputs, "G1_"));
                oNames.AddRange(G2.OutputNames ?? makeStringNames(G2.nOutputs, "G2_"));
            }

            var A = combineDiag(G1.A, G2.A);
            var B = combineDiag(G1.B, G2.B);
            var C = combineDiag(G1.C, G2.C);
            var D = combineDiag(G1.D, G2.D);

            var G = new StateSpaceDiscrete(A, B, C, D, copyData: false)
            {
                InputNames = iNames?.ToArray(),
                OutputNames = oNames?.ToArray(),
                StateNames = sNames?.ToArray(),
                Ts = G1.Ts,
            };

            return G;

            // todo: link G2's output to G1 input
            throw new NotImplementedException();
        }

        #endregion


        public IEnumerable<double[]> CalcResponse(IEnumerable<double[]> uVecs)
        {
            var uu = uVecs.Select(u => new DenseVector(u));
            return this.CalcResponse(uVecs);
        }

        
        public IEnumerable<Vector<double>> CalcResponse(IEnumerable<Vector<double>> uVecs)
        {
            var nSteps = uVecs.Count();
            // count should be complexity O(1) if a ICollection like input is given otherwise this statement will make the code somewhat slower
            var yVecs = new Vector<double>[nSteps];

            int i = 0;
            foreach (Vector<double> u in uVecs)
            {
                if (u == null)
                {
                    throw new ArgumentNullException("u vector at element " + i + " is null");
                }
                if (u.Count != nInputs)
                {
                    throw new ArgumentOutOfRangeException("the input vectors dimension did not align with the number of inputs expected (exp: {nInputs}, got: {u.Length})");
                }
                // calc output
                yVecs[i] = C * x_k1 + D * u;

                // calc updated states
                x_k1 = A * x_k1 + B * u;
                i++;
            }

            return yVecs;
        }

        public bool IsStable()
        {
            throw new NotImplementedException();
        }

        public IEnumerable<Vector<double>> Impulse(int nSteps)
        {
            throw new NotImplementedException();
        }

        public IEnumerable<Vector<double>> Impulse()
        {
            throw new NotImplementedException();
        }

        public IEnumerable<Complex[]> Bode(int nPoints = 100)
        {
            throw new NotImplementedException();
        }

        public IEnumerable<Complex[]> Bode(int nPoints, out double[] omega_vec)
        {
            throw new NotImplementedException();
        }

        public IEnumerable<Complex[]> Bode(double[] omega_vec)
        {
            throw new NotImplementedException();
        }

        public IEnumerable<Complex[]> GetPoles()
        {
            throw new NotImplementedException();
        }

        public IEnumerable<Complex[]> GetZeros()
        {
            throw new NotImplementedException();
        }

    }
}
