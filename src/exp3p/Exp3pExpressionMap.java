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

import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author tkonopka
 */
class Exp3pExpressionMap {

    final HashMap<String, ExpressionInfo> expinfo = new HashMap<String, ExpressionInfo>();
    final Exp3pAnnotationMap annomap;
    private final int numsamples;

    /**
     * Constructor remembers objects passed in constructor and sets up data
     * structures for recording expression values for all entries in the
     * annotation map and for all samples
     *
     * @param annomap
     * @param numsamples
     */
    public Exp3pExpressionMap(Exp3pAnnotationMap annomap, int numsamples) {
        this.annomap = annomap;
        this.numsamples = numsamples;

        // create expression storage space for each region defined in the annotation map
        ArrayList<String> temp = annomap.regionnames;               
        for (int i = 0; i < temp.size(); i++) {            
            expinfo.put(temp.get(i), new ExpressionInfo(numsamples));
        }
    }
    
    /**
     * allows to retrieve objects from the hashmap
     * 
     * @param name
     * @return 
     */
    public ExpressionInfo get(String name) {
        return expinfo.get(name);
    }
    
    public int getNumSamples() {
        return numsamples;
    }
}
