/*
 * Copyright 2013 Tomasz Konopka.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package exp3p;

import java.text.DecimalFormat;

/**
 * Container for expression information on a transcript.
 *
 * Designed to store information on multiple samples in one big table.
 *
 *
 * @author tkonopka
 */
public class ExpressionInfo {

    // info will holds values representative of expression
    // array is structed as [numsamples][4]
    // the 4 items as second dimension are total counts and counts weighted in various ways
    private final double[][] info;
    private final int numsamples;
    // indeces for the second dimension of the info array
    public static int TOT = 0;
    public static int WEIGHTEDMQ = 1;
    public static int WEIGHTEDTR = 2;
    public static int WEIGHTEDMQTR = 3;
    private final DecimalFormat scoreformat = new DecimalFormat("0.00");

    /**
     * initializes a matrix of doubles. By default all values will be 0.0
     *
     * @param numsamples
     */
    public ExpressionInfo(int numsamples) {
        this.numsamples = numsamples;
        info = new double[numsamples][4];
    }

    /**
     * This increments values associated with a sample.
     *
     * @param sample
     *
     * Must be an integer consistent with the number of samples declared in
     * constructor. Function does not check this.
     *
     * @param code
     *
     * Must be an integer from 0 to 3, see the constants declared in header.
     *
     * @param value
     *
     * The value to add to the specified sample and code.
     *
     */
    public void add(int sample, int code, double value) {
        info[sample][code] += value;
    }

    /**
     *
     * @return
     *
     * reports all the doubles in one row, tab separated (Has initial tab)
     *
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < numsamples; i++) {
            for (int j = 0; j < 4; j++) {
                sb.append("\t").append(scoreformat.format(info[i][j]));
            }
        }
        return sb.toString();
    }

        
    public double getInfoAt(int sample, int code) {
        return info[sample][code];
    }
}
