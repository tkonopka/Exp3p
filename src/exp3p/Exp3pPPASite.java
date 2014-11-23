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
 * holder class for a genomic position plus annotations 
 * 
 * @author tkonopka
 */
class Exp3pPPASite extends GenomicPosition {
    
    private String label="NA";
    private char strand='.';
    
    public Exp3pPPASite(String chr, int position) {
        super(chr, position);
    }

    public Exp3pPPASite(GenomicPosition gpos) {
        super(gpos.getChr(), gpos.getPosition());
    }
    
    public String getLabel() {
        return label;
    }

    public void setLabel(String label) {
        this.label = label;
    }

    public char getStrand() {
        return strand;
    }

    public void setStrand(char strand) {
        this.strand = strand;
    }    
    
    /**
     * 
     * @return 
     * 
     * a representation of the site suitable for writing in a bed file. 
     * 
     * Note: the site coordinates are shown in 0-based system in bed format,
     * although it is stored in 1-based system internally.
     * 
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(128);
        sb.append(getChr()).append("\t").append(this.getPosition()-1).append("\t").append(this.getPosition());
        sb.append("\t").append(strand).append(label);
        return sb.toString();
    }
}
