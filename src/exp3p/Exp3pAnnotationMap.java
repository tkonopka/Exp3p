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

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import jsequtils.file.BufferedReaderMaker;

/**
 * Object holds
 *
 *
 * @author tkonopka
 */
class Exp3pAnnotationMap {

    // annomap and fullannomap will hold regions of interest, e.g. genes and transcripts
    // the hashmaps will be accessed by chromosome names?
    final private HashMap<String, IntervalWithValueList> txsmap =
            new HashMap<String, IntervalWithValueList>(128);
    final private HashMap<String, IntervalWithValueList> fullannomap =
            new HashMap<String, IntervalWithValueList>(128);
    // mappings from transcript names to gene names and vice versa
    final HashMap<String, String> tx2gene = new HashMap<String, String>(128);
    final HashMap<String, StringBuilder> gene2tx = new HashMap<String, StringBuilder>(128);
    // regionnames will store the names of all transcripts, in the order in which they were defined
    final ArrayList<String> regionnames = new ArrayList<String>(1024);
    final ArrayList<String> genenames = new ArrayList<String>(1024);
    // txlengths will store, for a transcript name, the length of the region where counting
    // will be performed
    final private HashMap<String, Integer> txlength = new HashMap<String, Integer>();
    // fultxlength will hold the entire length of the transcript (as opposed only the 3prime length)
    final private HashMap<String, Integer> fulltxlength = new HashMap<String, Integer>();

    // Enum will enable to distinguish how information from an input file is going to be parsed in the
    // constructor
    static enum AnnoFileTypeEnum {

        REFGENE, BED
    }

    /**
     *
     * @param loadfile
     * @param filetype
     * @param numsamples
     * @param prime3len
     *
     * set this when loading REFGENE type input files
     *
     * @throws IOException
     */
    public Exp3pAnnotationMap(File loadfile, AnnoFileTypeEnum filetype, int prime3len, int extendby) throws IOException {
        switch (filetype) {
            case REFGENE:
                loadAnno(loadfile, prime3len);
                break;
            case BED:
                loadBed(loadfile, extendby);
                break;
            default:
                break;
        }
        // before exiting, make sure that intervals on all chromosomes are sorted
        sortAnnoComponents();
    }

    /**
     * function is called at the end of the constructor to make sure that the
     * intervals defined in the annotation objects are sorted by position
     *
     */
    private void sortAnnoComponents() {
        for (Map.Entry<String, IntervalWithValueList> entry : txsmap.entrySet()) {
            entry.getValue().sort();
        }
        // also in the full map
        for (Map.Entry<String, IntervalWithValueList> entry : fullannomap.entrySet()) {
            entry.getValue().sort();
        }
    }

    /**
     * Collects information about genes and adds them into annomap and
     * fullannomap. When finished, both the HashMaps should have intervals that
     * are properly sorted.
     *
     * @return
     * @throws IOException
     */
    private void loadAnno(File annofile, int prime3len) throws IOException {
        BufferedReader br = BufferedReaderMaker.makeBufferedReader(annofile);

        ArrayList<String> tempregionnames = new ArrayList<String>();
        ArrayList<String> tempgenenames = new ArrayList<String>();

        String s;
        while ((s = br.readLine()) != null) {
            if (!s.startsWith("#")) {
                String[] ss = s.split("\t");
                // call the transcript using the chromosome strand+official refseq name
                String transcript = ss[3] + ss[1];
                String chr = ss[2];

                IntervalWithValueList chromIntervals;

                // get a list of intervals for the three prime length            
                chromIntervals = txsmap.get(chr);
                if (chromIntervals == null) {
                    chromIntervals = new IntervalWithValueList();
                }
                int tx3primelen = get3pIntervalsFromAnno(transcript, chromIntervals, ss, prime3len);
                txsmap.put(chr, chromIntervals);
                txlength.put(transcript, tx3primelen);

                // get a list of intervals for the full transcript            
                chromIntervals = fullannomap.get(chr);
                if (chromIntervals == null) {
                    chromIntervals = new IntervalWithValueList();
                }
                int txfulllen = get3pIntervalsFromAnno(transcript, chromIntervals, ss, Integer.MAX_VALUE);
                fullannomap.put(chr, chromIntervals);
                fulltxlength.put(transcript, txfulllen);

                String genename = ss[12];

                tempregionnames.add(transcript);
                tempgenenames.add(genename);

                updateTranscriptGeneMaps(transcript, genename);
            }
        }

        // make sure that unique region names added in the loop are remembered
        setRegionNames(regionnames, tempregionnames);
        setRegionNames(genenames, tempgenenames);

        br.close();
    }

    /**
     * augments the tx2gene and gene2tx hashmaps with a new pair of
     * transcript/gene names
     *
     * @param transcript
     * @param genename
     */
    private void updateTranscriptGeneMaps(String transcript, String genename) {
        // transcript -> gene map update is easy because it is one to one
        tx2gene.put(transcript, genename);

        // gene -> transcript map update is more complicated because it is one to many
        if (gene2tx.containsKey(genename)) {
            StringBuilder nowtx = gene2tx.get(genename);
            nowtx.append(";").append(transcript);
            gene2tx.put(genename, nowtx);
        } else {
            StringBuilder nowtx = new StringBuilder();
            nowtx.append(transcript);
            gene2tx.put(genename, nowtx);
        }
    }

    /**
     * Collects information about genes and adds them into annomap and
     * fullannomap. When finished, both the HashMaps should have intervals that
     * are properly sorted.
     *
     * When reading the bedfile, intervals are converted into 1-based
     * coordinates.
     *
     * @return
     * @throws IOException
     */
    private void loadBed(File bedfile, int extendby) throws IOException {
        BufferedReader br = BufferedReaderMaker.makeBufferedReader(bedfile);

        ArrayList<String> tempregionnames = new ArrayList<String>();        
        
        String s;
        while ((s = br.readLine()) != null) {
            if (!s.startsWith("#")) {
                String[] ss = s.split("\t");
                String chr = ss[0];
                int start = 1 + Integer.parseInt(ss[1]);
                int end = Integer.parseInt(ss[2]);
                String id;
                if (ss.length > 3) {
                    if (ss[3].length() > 0 && (ss[3].charAt(0) == '+' || ss[3].charAt(0) == '-' || ss[3].charAt(0) == '.')) {
                        id = ss[3];
                    } else {
                        id = "." + ss[3];
                    }
                } else {
                    id = "." + ss[0] + ":" + ss[1] + "-" + ss[2];
                }

                IntervalWithValueList chromIntervals;

                // get a list of intervals also
                chromIntervals = txsmap.get(chr);
                if (chromIntervals == null) {
                    chromIntervals = new IntervalWithValueList();
                }

                // make an new interval that is slightly larger than officiallyl declared
                chromIntervals.add(start - extendby, end + extendby, id);
                txsmap.put(chr, chromIntervals);
                updateLength(txlength, id, end - start); // this used to have 1+end-start, was that right?

                // get a list of intervals for the full transcript            
                chromIntervals = fullannomap.get(chr);
                if (chromIntervals == null) {
                    chromIntervals = new IntervalWithValueList();
                }

                chromIntervals.add(start - extendby, end + extendby, id);
                fullannomap.put(chr, chromIntervals);
                updateLength(fulltxlength, id, end - start); // this used to have 1+end-start, was that right?

                // create space in memory where expression values will be recorded
                tempregionnames.add(id);
            }
        }

        // make sure that unique region names added in the loop are remembered
        setRegionNames(regionnames, tempregionnames);

        br.close();
    }

    /**
     * Adds intervals defined by annotation array (split from a refSeq record)
     * into the list aIntervals. Intervals of only up to maximym length
     * specified.
     *
     * Formally returns nothing, but the aIntervals object will be modified.
     *
     * @param aIntervals
     * @param astring
     * @param maxlen
     */
    private int get3pIntervalsFromAnno(String transcript, IntervalWithValueList aIntervals, String[] astring, int maxlen) {

        //String chr = astring[2];
        String strand = astring[3];
        String[] start = astring[9].split(",");
        String[] end = astring[10].split(",");

        // convert the start and end coordinates into true integers
        int[] startint = new int[start.length];
        int[] endint = new int[end.length];
        for (int i = 0; i < startint.length; i++) {
            startint[i] = 1 + Integer.parseInt(start[i]);
            endint[i] = Integer.parseInt(end[i]);
        }

        // add portions of exons into the chromIntervals data structure
        int nowlen = 0;
        int nowexon = 0;
        if (strand.equals("-")) {
            // add from the "first" (i.e. left-most) exon                
            while (nowlen < maxlen && nowexon < startint.length) {
                int exonlen = endint[nowexon] - startint[nowexon] + 1;
                // decide whether add a portion of the exon or the whole thing
                if (nowlen + exonlen > maxlen) {
                    aIntervals.add(startint[nowexon], startint[nowexon] + (maxlen - nowlen) - 1, transcript);
                    nowlen = maxlen;
                } else {
                    // add the whole thing
                    aIntervals.add(startint[nowexon], endint[nowexon], transcript);
                    nowlen += (endint[nowexon] - startint[nowexon]) + 1;
                }
                nowexon++;
            }
        } else {
            // on positive strand, add regions from the far right and move leftward
            nowexon = start.length - 1;
            while (nowexon >= 0 && nowlen < maxlen) {
                int exonlen = endint[nowexon] - startint[nowexon] + 1;
                if (nowlen + exonlen > maxlen) {
                    // add only a small part
                    aIntervals.add(endint[nowexon] - (maxlen - nowlen) + 1, endint[nowexon], transcript);
                    nowlen = maxlen;
                } else {
                    // add the whole thing
                    aIntervals.add(startint[nowexon], endint[nowexon], transcript);
                    nowlen += (endint[nowexon] - startint[nowexon]) + 1;
                }
                nowexon--;
            }
        }

        return nowlen;
    }

    private void setRegionNames(ArrayList<String> globnames, ArrayList<String> tempnames) {
        globnames.clear();
        Collections.sort(tempnames);
        if (tempnames.isEmpty()) {
            return;
        }
        globnames.add(tempnames.get(0));
        for (int i = 1; i < tempnames.size(); i++) {
            String nowname = tempnames.get(i);
            if (!nowname.equals(tempnames.get(i - 1))) {
                globnames.add(nowname);
            }
        }
    }

    /**
     * Inserts an id-length association into a hashmap. If an id is already
     * present in the map (i.e. there is a length association already), the
     * function stores the sum of the old and new length value.
     *
     * @param stringint
     * @param id
     * @param length
     */
    private void updateLength(HashMap<String, Integer> stringint, String id, int length) {
        if (stringint.containsKey(id)) {
            int oldint = stringint.get(id);
            stringint.put(id, length + oldint);
        } else {
            stringint.put(id, length);
        }
    }

    /**
     *
     * @param transcript
     * @return
     *
     * if ExpressionMap was loaded with an annotation table, this should return
     * the name of a gene associated with a transcript name. If the annotation
     * was not loaded, or the transcript name is incorrect, this will return
     * null
     *
     */
    public String getGeneFromTx(String transcript) {
        return tx2gene.get(transcript);
    }

    /**
     * 
     * @param genename
     * @return 
     * 
     * this should return an array with names of all transcripts associated with a gene name
     */
    public String[] getTxFromGene(String genename) {
        String txs = gene2tx.get(genename).toString();
        String[] tokens = txs.split(";");
        return tokens;
    }
    
    /**
     *
     * This is a defensive copy of an internal ArrayList.
     *
     * @return
     *
     * an array with all the names declared in the annotation object.
     *
     */
    public String[] getAllTxNames() {
        int alen = regionnames.size();
        String[] ans = new String[alen];
        for (int i = 0; i < alen; i++) {
            ans[i] = regionnames.get(i);
        }
        return ans;
    }

    /**
     * 
     * @return 
     * 
     * an array with all the names of genes declared in the annotation object
     * (loaded in loadAnno)
     * 
     */
    public String[] getAllGeneNames() {
        int alen = genenames.size();
        String[] ans = new String[alen];
        for (int i = 0; i < alen; i++) {
            ans[i] = genenames.get(i);
        }
        return ans;
    }
    
    /**
     *
     * @param chr
     * @return
     *
     * true if the underlying annotation map has information about a chromosome
     *
     */
    public boolean hasChr(String chr) {
        return txsmap.containsKey(chr);
    }

    /**
     *
     * @param chr
     * @return
     *
     * the intervals defined in the annotation files for a chromosome.
     *
     */
    public IntervalWithValueList getIntervals(String chr) {
        if (!hasChr(chr)) {
            return new IntervalWithValueList();
        }
        return txsmap.get(chr);
    }

    /**
     *
     * @param chr
     * @param position
     * @return
     *
     * a list of annotation labels associated with a genomic coordinate
     */
    public ArrayList<String> getValues(String chr, int position) {
        if (!hasChr(chr)) {
            return new ArrayList<String>(1);
        }
        return txsmap.get(chr).getValues(position);
    }

    /**
     *
     * Similar to getValues, but fetches labels from a different annotation
     * object (one with regions defined verbatim from the annotation file, not
     * from 3p ends)
     *
     * @param chr
     * @param position
     * @return
     */
    @Deprecated
    public ArrayList<String> getValuesFull(String chr, int position) {
        if (!hasChr(chr)) {
            return new ArrayList<String>(1);
        }
        return fullannomap.get(chr).getValues(position);
    }

    public int getTxLength(String name) {
        return txlength.get(name);
    }

    public int getFullTxLength(String name) {
        return fulltxlength.get(name);
    }
}
