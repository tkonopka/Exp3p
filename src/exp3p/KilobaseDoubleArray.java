/*
 * Copyright 2014 Tomasz Konopka.
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

/**
 * Similar to KilobaseIntArray, but stores doubles instead of ints.
 *
 *
 *
 * @author tkonopka
 */
public class KilobaseDoubleArray {

    private final double[][] values;
    private final int totlength;
    private final int numchunks;

    public KilobaseDoubleArray(int length) {
        // store the total length
        this.totlength = length;

        // assign an array of chunks
        int temp = length / 1024;
        if (length % 1024 > 0) {
            temp++;
        }
        numchunks = temp;
        values = new double[numchunks][];

        // set the intial values of the chunks to null, they will be initialized later when needed
        for (int i = 0; i < numchunks; i++) {
            values[i] = null;
        }
    }

    /**
     * Add one to the current value stored at an index
     *
     * @param index
     */
    public void increment(int index, double value) {
        // find the id of chunk
        int nowchunk = index / 1024;
        // the small index is equal to the last 10 bits of the index, which can be got by 
        // bitwise AND with a mask 000..001111111111.
        int smallindex = index & (1023);
        // set the new value at the position
        if (values[nowchunk] == null) {
            values[nowchunk] = new double[1024];
            values[nowchunk][smallindex] = value;
        } else {
            values[nowchunk][smallindex] += value;
        }
    }

    /**
     * Set a value at chosen locus of the array
     *
     * @param index
     * @param value
     */
    public void set(int index, double value) {
        // find the id of chunk
        int nowchunk = index / 1024;
        // the small index is equal to the last 10 bits of the index, which can be got by 
        // bitwise AND with a mask 000..001111111111.
        int smallindex = index & (1023);
        // set the new value at the position
        if (values[nowchunk] == null) {
            values[nowchunk] = new double[1024];
            values[nowchunk][smallindex] = value;
        } else {
            values[nowchunk][smallindex] = value;
        }
    }

    /**
     * obtain the value stored at a particular index of the array
     *
     * @param index
     * @return
     */
    public double get(int index) {
        // find the id of the chunk
        int nowchunk = index / 1024;
        if (values[nowchunk] == null) {
            return 0;
        } else {
            int smallindex = index - (nowchunk * 1024);
            return values[nowchunk][smallindex];
        }
    }

    /**
     *
     * @return
     *
     * the number of chunks that have been assigned in memory. This is a proxy
     * for the memory usage of the object.
     *
     *
     */
    public int getNumChunks() {
        int count = 0;
        for (int i = 0; i < numchunks; i++) {
            if (values[i] != null) {
                count++;
            }
        }
        return count;
    }

    /**
     *
     * @return
     *
     * the value originally used in the constructor, i.e. the number of elements
     * that can be stored in the array (actual number of elements may actually
     * be larger because chunks are always allocated by kb)
     */
    public int size() {
        return totlength;
    }
}
