package org.orangepalantir.leastsquares;

/** This code, by Matthew B. Smith, is available under an MIT licence
 * see https://github.com/odinsbane/least-squares-in-java/
 * 
 * We use this software here as part of a larger software, as permitted by the MIT license
 * The original licence header reads:
 * The MIT License (MIT)

Copyright (c) 2014 odinsbane

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * 
 */

public interface Function{
    /**
     *      Returns the functions evaluated at the specific parameter set
     * @return needs to evaluate the function
     * */
    public double evaluate(double[] values, double[] parameters);
    public int getNParameters();
    public int getNInputs();
    
}

