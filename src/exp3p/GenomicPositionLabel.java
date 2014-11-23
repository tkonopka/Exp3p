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
 * container class with chromosime, position, and a label (e.g. a gene name)
 * 
 * @author tkonopka
 */
public class GenomicPositionLabel extends GenomicPosition {

    String label;

    public GenomicPositionLabel(String chr, int position, String label) {
        super(chr, position);
        this.label = label;
    }

    public String getLabel() {
        return label;
    }
}
