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
 * 
 */
package exp3p;

/**
 *
 * Gateway to the Exp3p package/software suite. Tools are meant for analysis of RNA-seq data, 
 * particularly 3-prime seq with T-fill. Tools include evaluating expression from bam files, calling
 * 3-prime transcript end sites from bam files.
 * 
 * 
 * 
 * 
 * @author tkonopka
 */
public class Exp3p {

    final static String version = "0.1.0";

    public static String getVersion() {
        return version;
    }

    public static void printHelp() {
        System.out.println("Exp3p evaluates expression on genes, transcripts, or regions");
        System.out.println("using aligned read counts. It can also handle data from RNA-seq");
        System.out.println("protocols that measure reads only near transcript poly-A tails.");
        System.out.println();
        System.out.println("Usage: java -jar Exp3p.jar TYPE [options]");
        System.out.println();
        System.out.println("  eval          - evaluate coverage on transcripts or regions");
        System.out.println("  labelArich    - in bam file, label reads in AT rich regions (for Tfill)");
        System.out.println("  callPPA       - identify pre-poly-A sites");
        System.out.println("  richregions   - find regions rich in one nucleotide (e.g. A-rich regions)");
        System.out.println("  version       - display program version");
        System.out.println();
        System.out.println("\nAuthor: Tomasz Konopka (tkonopka@ludwig.ox.ac.uk)\n");
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        if (args == null || args.length == 0) {
            printHelp();
            return;
        }

        // copy all but one of the argument into a different arguments array
        String[] newargs = new String[args.length - 1];
        for (int i = 0; i < newargs.length; i++) {
            newargs[i] = args[i + 1];
        }

        // now call one of the tools in this toolkit
        if (args[0].equalsIgnoreCase("eval")) {
            new Exp3pEval(newargs).run();
        } else if (args[0].equalsIgnoreCase("labelArich")) {
            new Exp3pLabelArich(newargs).run();       
        } else if (args[0].equalsIgnoreCase("callPPA")) {
            new Exp3pCallPPA(newargs).run();        
        } else if (args[0].equalsIgnoreCase("version")) {
            System.out.println("Eval3p v" + version);
        } else if (args[0].equalsIgnoreCase("richregions")) {
            new RichRegions(newargs).run();        
        } else if (args[0].equalsIgnoreCase("help")) {
            printHelp();
        } else {
            System.out.println("Unrecognized command: " + args[0]);
        }

    }
}
