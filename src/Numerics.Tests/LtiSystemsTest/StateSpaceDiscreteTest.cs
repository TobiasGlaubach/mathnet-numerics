// <copyright file="AkimaSplineTest.cs" company="Math.NET">
// Math.NET Numerics, part of the Math.NET Project
// http://numerics.mathdotnet.com
// http://github.com/mathnet/mathnet-numerics
//
// Copyright (c) 2009-2016 Math.NET
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
// </copyright>

using MathNet.Numerics.SystemModels.LtiSystems;
using NUnit.Framework;
using System;

namespace MathNet.Numerics.UnitTests.LtiSystemsTest
{
    [TestFixture, Category("LTI")]
    public class StateSpaceDiscreteTest
    {
        [Test]
        public void TestConstructWrongInput()
        {
            Assert.Fail("Not implemented");
            // todo write test
        }

        [Test]
        public void TestLinkWrongInput()
        {
            Assert.Throws<ArgumentException>(() =>
            {
                var G1 = new StateSpaceDiscrete(2) { Ts = 1 };
                var G2 = new StateSpaceDiscrete(2) { Ts = 2 };
                
                var a = StateSpaceDiscrete.Link(G1, G2);
            }, "this should fail, since the sampling times of both models are different");
            Assert.Throws<NotSupportedException>(() =>
            {
                var G1 = new StateSpaceDiscrete(2) { Ts = 1, OutputNames = new string[] { "n1", "n2" } };
                var G2 = new StateSpaceDiscrete(2) { Ts = 1 };
                var a = StateSpaceDiscrete.Link(G1, G2);
            }, "this should fail, since the in and outut names of both models do not match");
            Assert.Throws<NotSupportedException>(() =>
            {
                var G1 = new StateSpaceDiscrete(2) {};
                var G2 = new StateSpaceDiscrete(2) {InputNames = new string[] { "n1", "n2" } };
                var a = StateSpaceDiscrete.Link(G1, G2);
            }, "this should fail, since the in and outut names of both models do not match");
            Assert.Throws<NotSupportedException>(() =>
            {
                var G1 = new StateSpaceDiscrete(2) {OutputNames = new string[] { "n1", "n2" } };
                var G2 = new StateSpaceDiscrete(2) {OutputNames = new string[] { "n3", "n4" } };
                var a = StateSpaceDiscrete.Link(G1, G2);
            }, "this should fail, since the in and outut names of both models do not match");

            Assert.Throws<NotSupportedException>(() =>
            {
                var G1 = new StateSpaceDiscrete(2) { OutputNames = new string[] { "n1", "n2" } };
                var G2 = new StateSpaceDiscrete(3) { OutputNames = new string[] { "n3", "n4" } };
                var a = StateSpaceDiscrete.Link(G1, G2);

            }, "this should fail, since the Dimensions of both models do not align");
        }

    }

}
