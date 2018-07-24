using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace MathNet.Numerics.LtiSystems
{
    /// <summary>
    /// a class for LTI state space models in the form
    /// x_k1 = A * x_k  + B * u
    /// y_k1 = C * x_k1 + D * u
    /// </summary>
    public class StateSpaceDiscrete
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
        public double Ts { get; private set; }

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


        /// <summary>
        /// 
        /// </summary>
        /// <param name="nStates"></param>
        /// <param name="nInputs"></param>
        /// <param name="nOutputs"></param>
        public StateSpaceDiscrete(int nStates, int nInputs = 1, int nOutputs = 1)
        {
            _A = DenseMatrix.CreateDiagonal(nStates, nStates, 1.0);
            _B = DenseMatrix.CreateDiagonal(nStates, nInputs, 1.0);
            _C = DenseMatrix.CreateDiagonal(nOutputs, nInputs, 1.0);
            _D = new DenseMatrix(nOutputs, nInputs);
            checkDimensions(A, B, C, D, A.ColumnCount, B.ColumnCount, C.ColumnCount);
        }

        public StateSpaceDiscrete(Matrix<double> A, Matrix<double> B, Matrix<double> C, Matrix<double> D)
        {
            checkDimensions(A, B, C, D, A.ColumnCount, B.ColumnCount, C.ColumnCount);
            _A = A;
            _B = B;
            _C = C;
            _D = D;
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

        public static StateSpaceDiscrete operator *(StateSpaceDiscrete G1, StateSpaceDiscrete G2)
        {
            // todo: link G2's output to G1 input
            throw new NotImplementedException();
        }

        
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


        #endregion

    }
}
