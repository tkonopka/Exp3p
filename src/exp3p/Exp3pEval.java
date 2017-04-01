/*
 * Copyright 2013-2017 Tomasz Konopka.
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

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import jsequtils.file.OutputStreamMaker;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

/**
 * Exp3pEval is the main tool in the Exp3p package. It is meant to scan a bam
 * file and collect readcount information on regions of interest (defined
 * through bed or refSeq-style annotation files).
 *
 *
 *
 * @author tkonopka
 */
public class Exp3pEval implements Runnable {

    Exp3pStrandedEnum stranded = Exp3pStrandedEnum.NONE;
    // variables with values determined from the command line
    private int prime3len = 150;
    private double rescaleIntervalFactor = 1.0;
    private String command;
    private ArrayList<File> bamfiles;
    private ArrayList<String> labels;
    // when the bam files will be scanned, the number of mapped reads (in millions) 
    // will be stored in mappedreadsM
    private ArrayList<Double> mappedreadsM = new ArrayList<Double>();
    private File annofile = null, bedfile = null;
    //private File ppafile = null;
    private int bedextend = 4;
    //private File genomefile;
    private String output = "stdout";
    private String outputgene = null;
    // tfillnorm determines how expression values will be normalized
    // set tfillnorm to -1.0 to normalize by length as in FPKM or to a finite value
    // to use a constant for the length in FPKM
    private double tfillnorm = -1.0;
    // ignoreAR - set to true to ignore reads that labeled by the AR tag (false positives due to AT rich region)
    private boolean filterAR = false;
    // isok determines if class is initialized correctly
    private boolean isok = false;

    public void printEvalHelp() {
        System.out.println("Exp3p eval: evaluates expression on genes, transcripts, or regions using alignment counts");
        System.out.println();
        System.out.println("Usage: java -jar Exp3p.jar eval [options]");
        System.out.println();
        System.out.println("  --bam <File>          - alignment file to process ");
        System.out.println("                          (can specify multiple times)");
        System.out.println("  --bed <File>          - regions file (alternative to --anno and --prime3len)");
        System.out.println("  --bedextend <int>     - extend bed intervals [default " + bedextend + "]");
        System.out.println("  --label <String>      - sample identifier");
        System.out.println("                          (must specify one label for each bam file)");
        System.out.println("  --anno <File>         - refGene-style annotation of transcripts");
        System.out.println("  --filterAR <boolean>  - filter out reads labeld by AR tag [default false]");
        System.out.println("  --prime3len <int>     - length of region on 3prime end of transcript");
        System.out.println("                          (default 150, set to 0 to consider full transcripts)");
        System.out.println("  --rescale <double>    - rescale the low-high expression estimate by a factor [default 1.0]");
        System.out.println("  --stranded <String>   - stranded protocol, use NONE (default), SAME, or REVERSE");
        System.out.println("  --output <String>     - output file will contain transcript-wise expression values");
        System.out.println("  --outputgene <String> - output file will contain gene-wise expression values");
        System.out.println("  --tfillnorm <double>  - normalization for tfill protocol [default -1 for FPKM normalization]");
        System.out.println();
        System.out.println("\nAuthor: Tomasz Konopka (tkonopka@ludwig.ox.ac.uk)\n");
    }

    private boolean parseEvalParameters(String[] args) {
        OptionParser prs = new OptionParser();

        // bam - input bam file
        prs.accepts("bam").withRequiredArg().ofType(File.class);
        // samplelabel - samplelabel that appear in vcf file
        prs.accepts("label").withRequiredArg().ofType(String.class);

        // anno - input annotation file
        prs.accepts("anno").withRequiredArg().ofType(File.class);
        // prime3len - length of 3prime end to look at
        prs.accepts("prime3len").withRequiredArg().ofType(Integer.class);

        // rescale - rescales low-high interval by a certain factor
        prs.accepts("rescale").withRequiredArg().ofType(Double.class);

        // stranded - stranded protocol
        prs.accepts("stranded").withRequiredArg().ofType(String.class);

        // alternative to the annotation file, also accept bed files
        prs.accepts("bed").withRequiredArg().ofType(File.class);
        prs.accepts("bedextend").withRequiredArg().ofType(Integer.class);

        // genome - genome annotation file
        prs.accepts("genome").withRequiredArg().ofType(File.class);

        // output - output file
        prs.accepts("output").withRequiredArg().ofType(String.class);
        prs.accepts("outputgene").withRequiredArg().ofType(String.class);

        // tfillnorm - normalization in lieu of transcript length
        prs.accepts("tfillnorm").withRequiredArg().ofType(Double.class);
        prs.accepts("filterAR").withRequiredArg().ofType(Boolean.class);

        // now use OptionSet to parse the command line
        OptionSet options;
        try {
            options = prs.parse(args);
        } catch (Exception ex) {
            System.out.println("Error parsing command line parameters\n" + ex.getMessage());
            return false;
        }

        // determine how much of the of the specified regions should be counted
        if (options.has("prime3len")) {
            prime3len = (int) (Integer) options.valueOf("prime3len");
            if (prime3len < 0) {
                System.out.println("Parameter prime3len must be non negative");
                return false;
            }
            // interpret prime3len=0 as request to use full transcripts as opposed to 3prime regions
            if (prime3len == 0) {
                prime3len = Integer.MAX_VALUE;
            }
        }

        // determine the regions that the user wants to look at
        if (options.has("anno") || options.has("bed") || options.has("ppa")) {
            int filecount = 0;
            if (options.has("anno")) {
                filecount++;
            }
            if (options.has("bed")) {
                filecount++;
            }
            if (filecount > 1) {
                System.out.println("Cannot use more than one of --bed --anno --ppa");
                return false;
            }
            if (options.has("anno")) {
                annofile = (File) options.valueOf("anno");
                if (!annofile.exists() || !annofile.canRead()) {
                    System.out.println("Cannot read annotation file: " + annofile.getAbsolutePath());
                    return false;
                }
            }
            if (options.has("bed")) {
                bedfile = (File) options.valueOf("bed");
                if (!bedfile.exists() || !bedfile.canRead()) {
                    System.out.println("Cannot read bed file: " + bedfile.getAbsolutePath());
                    return false;
                }
            }
        } else {
            System.out.println("missing required parameter bed, anno, or ppa ");
            return false;
        }

        // choose strandedness protocol
        if (options.has("stranded")) {
            String strandedString = (String) options.valueOf("stranded");
            if (strandedString.equalsIgnoreCase("NONE")) {
                stranded = Exp3pStrandedEnum.NONE;
            } else if (strandedString.equalsIgnoreCase("REVERSE")) {
                stranded = Exp3pStrandedEnum.REVERSE;
            } else if (strandedString.equalsIgnoreCase("SAME")) {
                stranded = Exp3pStrandedEnum.SAME;
            } else {
                System.out.println("Unrecognized value for parameter stranded: " + strandedString);
                return false;
            }
        }

        if (options.has("rescale")) {
            rescaleIntervalFactor = (Double) options.valueOf("rescale");
            if (rescaleIntervalFactor < 0) {
                System.out.println("Interval rescale factor must be > 0");
                return false;
            }
        }

        // check for inputs for alignment/labels
        if (options.has("bam")) {
            bamfiles = new ArrayList<File>((List<File>) options.valuesOf("bam"));
        } else {
            System.out.println("Missing required argument: bam");
            return false;
        }
        if (options.has("label")) {
            labels = new ArrayList<String>((List<String>) options.valuesOf("label"));
        } else {
            System.out.println("Missing required argument: label");
            return false;
        }

        // the number of input files and labels must match to uniquely associated         
        if (labels.size() != bamfiles.size()) {
            System.out.println("number of bam and label arguments do not match");
            return false;
        }

        // figure out where to write the output
        if (options.has("output")) {
            output = (String) options.valueOf("output");
        }
        if (options.has("outputgene")) {
            outputgene = (String) options.valueOf("outputgene");
        }

        if (options.has("tfillnorm")) {
            tfillnorm = (Double) options.valueOf("tfillnorm");
            if (tfillnorm == 0.0) {
                System.out.println("tfillnorm must be set to a non-zero number (use -1 to disable tfill normalization)");
                return false;
            }
        }
        if (options.has("filterAR")) {
            filterAR = (Boolean) options.valueOf("filterAR");
        }

        if (options.has("bedextend")) {
            bedextend = (Integer) options.valueOf("bedextend");
            if (bedextend < 0) {
                System.out.println("bedextend must be set >= 0");
                return false;
            }
        }

        return true;
    }

    
    /**
     * Constructor. Parses command line parameters only. The main computation 
     * work is performed in run().
     * 
     * 
     * @param args 
     * 
     * Command line parameters passed to this tool
     * 
     */
    public Exp3pEval(String[] args) {
        if (args == null || args.length == 0) {
            printEvalHelp();
            return;
        }

        if (parseEvalParameters(args)) {
            isok = true;
        }

        // record the command given into a variable
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < args.length; i++) {
            sb.append(" ").append(args[i]);
        }
        command = sb.toString();
    }

    /**
     * Execute the expression counting tool. Assumes the object is properly
     * initialized. 
     * 
     */
    @Override
    public void run() {
        if (!isok) {
            return;
        }

        // annomap and expmap will hold annotation and expression information, respectively
        Exp3pAnnotationMap annomap = null;
        Exp3pExpressionMap expmap = null;
        try {
            if (annofile != null) {
                annomap = new Exp3pAnnotationMap(annofile, Exp3pAnnotationMap.AnnoFileTypeEnum.REFGENE, prime3len, 0);
            } else if (bedfile != null) {
                annomap = new Exp3pAnnotationMap(bedfile, Exp3pAnnotationMap.AnnoFileTypeEnum.BED, prime3len, bedextend);
            }
            expmap = new Exp3pExpressionMap(annomap, bamfiles.size());
        } catch (Exception ex) {
            System.out.println("Exception while loading annotation: " + ex.getMessage());
            return;
        }

        // now loop through all the bam files and collect relevant expression info
        for (int i = 0; i < bamfiles.size(); i++) {
            // add an element to the mappedreadsM array. For now, this is set to 0.0
            // but next command processOneBam replaces this value with the true number
            mappedreadsM.add(0.0);
            processOneBam(annomap, expmap, bamfiles.get(i), i);
        }

        // print out all the recorded expression information                 
        try {
            if (annofile != null) {
                outputExpInfo(annomap, expmap);
                outputGeneExpInfo(annomap, expmap);
            } else if (bedfile != null) {
                outputExpInfo(annomap, expmap);
            }
        } catch (Exception ex) {
            System.out.println("Exception while outputing expression values: " + ex.getMessage());
        }

    }

    /**
     * This function parses one read, determines where the read is (is it on
     * region of interest?) and updates relevant expression counters.
     *
     * @param annomap
     * @param expmap
     * @param sampleindex
     * @param txs
     *
     * names of transcripts. The function assumes the transcript string will be
     * formatted with the strand as the first character, followed by the actual
     * transcript name, e.g. "-TP53" because TP53 is on the negative strand.
     *
     * @param record
     */
    private void updateExpInfo(Exp3pExpressionMap expmap,
            int sampleindex, ArrayList<String> txs, SAMRecord record, int valNH) {

        //System.out.println(record.getReferenceName()+":"+record.getAlignmentStart());
        if (txs.isEmpty()) {
            return;
        }

        int numtxs = txs.size();
        boolean[] countthis = new boolean[numtxs];

        // get number of transcripts that match strandedness criteria
        int numstranded = 0;
        if (stranded == Exp3pStrandedEnum.NONE) {
            for (int i = 0; i < numtxs; i++) {
                countthis[i] = true;
            }
        } else if (stranded == Exp3pStrandedEnum.SAME) {
            char goodstrand;
            if (record.getReadNegativeStrandFlag()) {
                goodstrand = '-';
            } else {
                goodstrand = '+';
            }
            for (int i = 0; i < numtxs; i++) {
                countthis[i] = (txs.get(i).charAt(0) == goodstrand) || (txs.get(i).charAt(0) == '.');
            }
        } else {
            char goodstrand;
            if (record.getReadNegativeStrandFlag()) {
                goodstrand = '+';
            } else {
                goodstrand = '-';
            }
            for (int i = 0; i < numtxs; i++) {
                countthis[i] = (txs.get(i).charAt(0) == goodstrand) || (txs.get(i).charAt(0) == '.');
            }
        }

        // count number of regions declared on plus and minus strands                            
        for (int i = 0; i < numtxs; i++) {
            if (countthis[i]) {
                numstranded++;
            }
        }

        // do not modify anything if none of the transcripts satisfy strandedness criteria
        if (numstranded == 0) {
            return;
        }

        for (int j = 0; j < numtxs; j++) {
            if (countthis[j]) {
                ExpressionInfo nowinfo = expmap.get(txs.get(j));
                nowinfo.add(sampleindex, ExpressionInfo.TOT, 0.5);
                nowinfo.add(sampleindex, ExpressionInfo.WEIGHTEDMQ, 0.5 / valNH);
                nowinfo.add(sampleindex, ExpressionInfo.WEIGHTEDTR, 0.5 / valNH);
                nowinfo.add(sampleindex, ExpressionInfo.WEIGHTEDMQTR, 0.5 / (valNH * numstranded));
            }
        }

    }

    /**
     * This is where a bam file is scanned and expression in regions of interest
     * are recorded into the "expinfo" object.
     *
     * @param annomap 
     * object holding regions of interest
     * @param expmap 
     * object holding expression estimates
     * @param infile     
     * File to process
     * @param sampleindex
     * index of samples (associated content of File to the expmap object)
     * 
     */
    private void processOneBam(Exp3pAnnotationMap annomap, Exp3pExpressionMap expmap, File infile, int sampleindex) {
        // set up input/output objects using SAM library
        SAMFileReader inputSam = new SAMFileReader(infile);
        inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        // count the number of aligned reads in the file
        double nummapped = 0.0;

        // keep track of current chromosome and position with the next 
        int lastchr = -1;
        int laststart = -1;
        int lastend = -1;
        ArrayList<String> starttxs = null;
        ArrayList<String> endtxs = null;

        for (final SAMRecord record : inputSam) {
            int nowchr = record.getReferenceIndex();
            int startpos = record.getAlignmentStart();
            int endpos = record.getAlignmentEnd();

            // figure out whether to use this read or not, based on its AR tag            
            boolean useThisRead = true && !record.getReadUnmappedFlag();
            if (filterAR) {
                Object ARtag = record.getAttribute("AR");
                if (ARtag != null) {
                    if ((Integer) ARtag == 1) {
                        useThisRead = false;
                    }
                }
            }

            if (useThisRead) {
                int valNH = 1;
                if (record.getAttribute("NH") != null) {
                    valNH = record.getIntegerAttribute("NH");
                }

                // record the start of the read
                // check if list of transcript at this locus should be reloaded            
                if (nowchr != lastchr || startpos != laststart || starttxs == null) {
                    String nowchrname = record.getReferenceName();
                    starttxs = annomap.getValues(nowchrname, startpos);
                }
                updateExpInfo(expmap, sampleindex, starttxs, record, valNH);

                // record the end of the read            
                if (nowchr != lastchr || endpos != lastend || endtxs == null) {
                    String nowchrname = record.getReferenceName();
                    endtxs = annomap.getValues(nowchrname, endpos);
                }
                updateExpInfo(expmap, sampleindex, endtxs, record, valNH);

                if (nowchr > -1) {
                    nummapped += 1.0 / valNH;
                }
            }

            lastchr = nowchr;
            laststart = startpos;
        }

        // at the end, close all the IO objects
        inputSam.close();

        // record the number of mapped reads into the main array                
        mappedreadsM.set(sampleindex, new Double(((double) nummapped) / 1e6));
    }

    /**
     *
     *
     * @return
     *
     * A string that contains multiple lines, all beginning with # (as a comment
     * on a table) The comment contains information about the samples
     *
     *
     *
     */
    private String getOutputHeader() {
        StringBuilder sb = new StringBuilder();
        SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd hh:mm:ss");


        sb.append("# exp3p eval ").append("\n");
        sb.append("# ").append(sdf.format(new Date())).append("\n");
        sb.append("#\n");
        sb.append("# command:\teval ").append(command).append("\n");
        sb.append("#\n");
        sb.append("# prime3len\t").append(prime3len).append("\n");
        if (bedfile == null) {
            sb.append("# annotation\t").append(annofile.getAbsolutePath()).append("\n");
        } else {
            sb.append("# bed\t").append(bedfile.getAbsolutePath()).append("\n");
        }
        sb.append("# rescaleIntervalFactor\t").append(rescaleIntervalFactor).append("\n");
        sb.append("#\n");
        int numsamples = bamfiles.size();
        for (int i = 0; i < numsamples; i++) {
            sb.append("# Sample\t").append(labels.get(i)).append("\t").append(bamfiles.get(i).getAbsolutePath()).append("\n");
        }
        sb.append("#\n");

        return sb.toString();
    }

    /**
     * Creates the final table with expression counts for all samples. Info in
     * the table comes from the object "expinfo" and goes into a stream defined
     * by "output"
     *
     * @param annomap
     * @param expmap
     * 
     * @throws FileNotFoundException
     * @throws IOException
     */
    private void outputExpInfo(Exp3pAnnotationMap annomap, Exp3pExpressionMap expmap) throws FileNotFoundException, IOException {

        StringBuilder sb;

        // create output stream
        OutputStream os = OutputStreamMaker.makeOutputStream(output);

        // output the header
        sb = new StringBuilder();
        sb.append(getOutputHeader());
        sb.append("Region.id\tStrand\tLength.full\tLength.3p");
        int numsamples = bamfiles.size();
        for (int i = 0; i < numsamples; i++) {
            String nowlabel = labels.get(i);
            sb.append("\t").append(nowlabel).append(".tot").
                    append("\t").append(nowlabel).append(".wmap").
                    append("\t").append(nowlabel).append(".wtxs").
                    append("\t").append(nowlabel).append(".wmaptxs").
                    append("\t").append(nowlabel).append(".EPM").
                    append("\t").append(nowlabel).append(".EPM.low").
                    append("\t").append(nowlabel).append(".EPM.high");
        }
        sb.append("\n");
        os.write(sb.toString().getBytes());

        // output information on all the transcripts
        DecimalFormat scoreformat1 = new DecimalFormat("0.000");
        DecimalFormat scoreformat2 = new DecimalFormat("0.###E0");
        sb = new StringBuilder();
        for (int i = 0; i < annomap.regionnames.size(); i++) {
            String nowname = annomap.regionnames.get(i);

            // get the length of the transcript region measured
            // this will be used to normalize the expression value by the length of transcript
            int nowlen = annomap.getTxLength(nowname);
            double nowlendbl = (double) nowlen;
            // but in tfill normalization, always renormalize by a fixed length, e.g. length of read
            if (tfillnorm > 0) {
                nowlendbl = tfillnorm;
            }

            sb.append(nowname.substring(1)).
                    append("\t").append(nowname.charAt(0)).
                    append("\t").append(annomap.getFullTxLength(nowname)).
                    append("\t").append(nowlen);

            ExpressionInfo ei = expmap.get(nowname);

            // write a row of information about each sample, each category of information
            for (int j = 0; j < numsamples; j++) {
                for (int k = 0; k < 4; k++) {
                    sb.append("\t").append(scoreformat1.format(ei.getInfoAt(j, k)));
                }
                // compute the EPM and low/high intervals
                double reads = ei.getInfoAt(j, 3);
                double readsLo = Math.max(0, reads - (rescaleIntervalFactor * Math.sqrt(reads + 0.5)));
                double readsHi = reads + (rescaleIntervalFactor * Math.sqrt(reads + 0.5));
                double totmapped = mappedreadsM.get(j);
                // adjust the totmapped by the length of the region
                totmapped = totmapped * nowlendbl / 1000.0;
                sb.append("\t").append(scoreformat2.format(reads / totmapped)).
                        append("\t").append(scoreformat2.format(readsLo / totmapped)).
                        append("\t").append(scoreformat2.format(readsHi / totmapped));
            }
            sb.append("\n");
            if (sb.length() > 100000) {
                os.write(sb.toString().getBytes());
                sb = new StringBuilder();
            }
        }
        // write leftover information
        if (sb.length() > 0) {
            os.write(sb.toString().getBytes());
        }

        // always clean up the streams at the end.
        if (os != System.out) {
            os.close();
        }
    }

    /**
     * Creates the final table with expression counts for all samples. Info in
     * the table comes from the object "expinfo" and goes into a stream defined
     * by "output"
     *
     * @param annomap
     * @param expmap
     * 
     * @throws FileNotFoundException
     * @throws IOException
     */
    private void outputGeneExpInfo(Exp3pAnnotationMap annomap, Exp3pExpressionMap expmap) throws FileNotFoundException, IOException {

        if (outputgene == null && annomap.genenames.size() > 0) {
            return;
        }

        StringBuilder sb;

        // create output stream
        OutputStream os = OutputStreamMaker.makeOutputStream(outputgene);

        // output the header
        sb = new StringBuilder();
        sb.append(getOutputHeader());
        sb.append("Region.id\tStrand");
        int numsamples = bamfiles.size();
        for (int i = 0; i < numsamples; i++) {
            String nowlabel = labels.get(i);            
            sb.append("\t").append(nowlabel).append(".count").
                    append("\t").append(nowlabel).append(".EPM").
                    append("\t").append(nowlabel).append(".EPM.low").
                    append("\t").append(nowlabel).append(".EPM.high");
        }
        sb.append("\n");
        os.write(sb.toString().getBytes());

        // output information on all the transcripts
        DecimalFormat scoreformat1 = new DecimalFormat("0.0");
        DecimalFormat scoreformat2 = new DecimalFormat("0.###E0");
        sb = new StringBuilder();
        for (int i = 0; i < annomap.genenames.size(); i++) {
            String nowgene = annomap.genenames.get(i);
            String[] alltxs = annomap.getTxFromGene(nowgene);

            // write the first few columns in output
            sb.append(nowgene).
                    append("\t").append(getConsensusStrand(alltxs));
            
            // write a row of information about each sample, each category of information
            for (int j = 0; j < numsamples; j++) {
                double reads = 0.0;
                double totmapped = mappedreadsM.get(j);

                double[] eest = new double[alltxs.length];
                double[] eestLo = new double[alltxs.length];
                double[] eestHi = new double[alltxs.length];

                double totSqDevLo = 0.0;
                double totSqDevHi = 0.0;
                double totEPM = 0.0;
                double totCount = 0.0;
                
                // collect readcounts on all transcripts
                for (int z = 0; z < alltxs.length; z++) {
                    String nowname = alltxs[z];
                    ExpressionInfo ei = expmap.get(nowname);
                    // get the length of the transcript region measured
                    // this will be used to normalize the expression value by the length of transcript
                    int nowlen = annomap.getTxLength(nowname);
                    double nowlendbl = (double) nowlen;
                    // but in tfill normalization, always renormalize by a fixed length, e.g. length of read
                    if (tfillnorm > 0) {
                        nowlendbl = tfillnorm;
                    }
                    // adjust the read counts   
                    reads = ei.getInfoAt(j, 3);
                    totmapped = totmapped * nowlendbl / 1000.0;
                    eest[z] = reads / totmapped;
                    eestLo[z] = Math.max(0, reads - (rescaleIntervalFactor * Math.sqrt(reads + 0.5))) / totmapped;
                    eestHi[z] = (reads + (rescaleIntervalFactor * Math.sqrt(reads + 0.5))) / totmapped;

                    // keep track of totals
                    totEPM += eest[z];
                    totSqDevLo += Math.pow(eestLo[z] - eest[z], 2);
                    totSqDevHi += Math.pow(eestHi[z] - eest[z], 2);
                    totCount += reads;
                }

                sb.append("\t").append(scoreformat1.format(totCount)).
                        append("\t").append(scoreformat2.format(totEPM)).
                        append("\t").append(scoreformat2.format(Math.max(0, totEPM - Math.sqrt(totSqDevLo)))).
                        append("\t").append(scoreformat2.format(totEPM + Math.sqrt(totSqDevHi)));
            }
            sb.append("\n");

            if (sb.length() > 100000) {
                os.write(sb.toString().getBytes());
                sb = new StringBuilder();
            }

        }
        // write leftover information
        if (sb.length() > 0) {
            os.write(sb.toString().getBytes());
        }

        // always clean up the streams at the end.
        if (os != System.out) {
            os.close();
        }
    }

    /**
     * Looks at the first characters of the strings and returns a consensus
     * strand.
     *
     *
     * @param txs
     * @return
     *
     * either '+' or '-' if all the transcripts give the same strand. Else
     * returns '.' if two transcripts disagree.
     *
     */
    private char getConsensusStrand(String[] txs) {
        char ans = txs[0].charAt(0);
        for (int i = 1; i < txs.length; i++) {
            if (ans != txs[i].charAt(0)) {
                return '.';
            }
        }
        return ans;
    }
}
