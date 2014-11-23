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

import jsequtils.genome.GenomicPosition;

/**
 * this is an object that holds a chromosome, start and end position.
 *
 * NOTE: since this is an extension of genomicposition, it has functions to set
 * "position". This class also allows to set a "start" and "end". The convention
 * is as follows. The start and end points are from 0-based coordinate system.
 * The position will be used in a 1-based coordinate system.
 *
 * E.g. a bed interval chr1 100 200 will return getStart() -> 100, getEnd() ->
 * 200, and getPosition() -> 101.
 *
 *
 *
 * @author tkonopka
 */
public class GenomicInterval extends GenomicPosition {

    int end;

    /**
     * constructor for new interval on one chromosome with start and end positions
     * 
     * @param chr
     * @param start
     * @param end 
     */
    public GenomicInterval(String chr, int start, int end) {
        super(chr, start + 1);
        this.end = end;
    }

    /**
     * copy constructor
     * 
     * @param gi 
     */
    public GenomicInterval(GenomicInterval gi) {
        super(gi.getChr(), gi.getPosition());
        this.end = gi.getEnd();
    }
    
    @Override
    public String toString() {
        return this.getChr() + "\t" + this.getStart() + "\t" + end;
    }

    @Override
    public String getChr() {
        return super.getChr();
    }

    @Override
    public void setChr(String chr) {
        super.setChr(chr);
    }

    public int getStart() {
        return super.getPosition() - 1;
    }

    public void setStart(int start) {
        super.setPosition(start + 1);
    }
          
    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public int getLength() {
        return end - getStart();
    }
    
    /**
     * Attempts to merge this interval with another one.
     * 
     * e.g. if this describes (bed) chr1 100 200
     * and the gi object describes (bed) chr1 150 250
     * then the merged object will be chr1 100 250 and returns true
     * 
     * @param gi
     * 
     * an interval to merge with
     * 
     * @return 
     * 
     * true if merge is succesful, false if not. 
     * 
     * When the return value is true, the object will be modified, i.e. made wider
     * to encompass both the old interval and the interval specified in the input.
     * 
     */
    public boolean mergeWith(GenomicInterval gi) {
        
        if (!getChr().equals(gi.getChr())) {
            return false;
        }
        
        int thisStart = getStart();
        int thisEnd = getEnd();
        int giStart = gi.getStart();
        int giEnd = gi.getEnd();
        if (giStart> thisEnd || thisStart > giEnd ) {
            return false;
        }
        
        System.out.println("merging: "+this.toString()+" with "+gi.toString());
        setStart(Math.min(thisStart, giStart));
        setEnd(Math.max(thisEnd, giEnd));
        
        System.out.println("to give: "+this.toString());
        return true;
    }
}
