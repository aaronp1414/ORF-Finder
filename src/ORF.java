import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class ORF {
    public static HashMap<String, String> acidMap5to3 = new HashMap<>();
    public static HashMap<String, String> acidMap3to5 = new HashMap<>();

    /*
        Initialize Static maps used for translation
        Idea to use HashMap for translation from:
        https://pdfs.semanticscholar.org/6aca/38fc20459010efdad3b535e71f43d8b17bc9.pdf
     */
    static {
        String[] AMINO_ACIDS = {
                "F", "F", "L", "L", "S",
                "S", "S", "S", "Y", "Y",
                "--STOP--", "--STOP--",
                "C", "C", "--STOP--",
                "W", "L", "L", "L", "L",
                "P", "P", "P", "P", "H",
                "H", "Q", "Q", "R", "R",
                "R", "R", "I", "I", "I", "M",
                "T", "T", "T", "T", "N", "N",
                "K", "K", "S", "S", "R", "R",
                "V", "V", "V", "V", "A", "A",
                "A", "A", "D", "D", "E", "E",
                "G", "G", "G", "G"
        };
        String[] CODONS5to3 = {
                "TTT", "TTC", "TTA",
                "TTG", "TCT", "TCC",
                "TCA", "TCG", "TAT",
                "TAC", "TAA", "TAG",
                "TGT", "TGC", "TGA",
                "TGG", "CTT", "CTC",
                "CTA", "CTG", "CCT",
                "CCC", "CCA", "CCG",
                "CAT", "CAC", "CAA",
                "CAG", "CGT", "CGC",
                "CGA", "CGG", "ATT",
                "ATC", "ATA", "ATG",
                "ACT", "ACC", "ACA",
                "ACG", "AAT", "AAC",
                "AAA", "AAG", "AGT",
                "AGC", "AGA", "AGG",
                "GTT", "GTC", "GTA",
                "GTG", "GCT", "GCC",
                "GCA", "GCG", "GAT",
                "GAC", "GAA", "GAG",
                "GGT", "GGC", "GGA",
                "GGG"
        };
        String[] CODONS3to5 = {
                "AAA", "AAG", "AAT",
                "AAC", "AGA", "AGG",
                "AGT", "AGC", "ATA",
                "ATG", "ATT", "ATC",
                "ACA", "ACG", "ACT",
                "ACC", "GAA", "GAG",
                "GAT", "GAC", "GGA",
                "GGG", "GGT", "GGC",
                "GTA", "GTG", "GTT",
                "GTC", "GCA", "GCG",
                "GCT", "GCC", "TAA",
                "TAG", "TAT", "TAC",
                "TGA", "TGG", "TGT",
                "TGC", "TTA", "TTG",
                "TTT", "TTC", "TCA",
                "TCG", "TCT", "TCC",
                "CAA", "CAG", "CAT",
                "CAC", "CGA", "CGG",
                "CGT", "CGC", "CTA",
                "CTG", "CTT", "CTC",
                "CCA", "CCG", "CCT",
                "CCC"
        };

        for (int i = 0; i < CODONS5to3.length; i++)
            acidMap5to3.put(CODONS5to3[i], AMINO_ACIDS[i]);
        for (int i = 0; i < CODONS3to5.length; i++)
            acidMap3to5.put(CODONS3to5[i], AMINO_ACIDS[i]);
    }

    private int startIndex;
    private int stopIndex;
    // sequence is always in 5' to 3' direction.
    private String sequence;
    private String acidSequence = null;
    // 0 meaning 5' to 3'. 1 meaning 3' to 5'
    private int direction;
    private ORF(int startIndex, int stopIndex, String sequence, int direction) {
        this.startIndex = startIndex;
        this.stopIndex = stopIndex;
        this.sequence = sequence;
        this.direction = direction;
    }

    /*
        Method to Find and Combine ORFs from given dna sequence in both 5' -> 3' and 3' -> 5' direction
         @param dnaSequence sequence to use
     */
    public static ArrayList<ORF> findORFs(StringBuilder dnaSequence) {
        ArrayList<ORF> allORFs = new ArrayList<>();
        //Find start and Stop codons in 5' to 3' Direction
        ArrayList<Integer> possibleStartCodonIndexes5to3 = findCodonIndexes(dnaSequence, "ATG");
        ArrayList<Integer> possibleStopCodonIndexes5to3 = findPossibleStopCodonIndexes5to3(dnaSequence);
        //Create the ORFS with the start and stop codons found
        ArrayList<ORF> orFs5to3 = createORFs5to3(dnaSequence, possibleStartCodonIndexes5to3, possibleStopCodonIndexes5to3);

        //Find start and Stop codons in 3' to 5' Direction
        ArrayList<Integer> possibleStartCodonIndexes3to5 = findCodonIndexes(dnaSequence, "CAT");
        ArrayList<Integer> possibleStopCodonIndexes3to5 = findPossibleStopCodonIndexes3to5(dnaSequence);
        //Create the ORFS with the start and stop codons found
        ArrayList<ORF> orFs3to5 = createORFs3to5(dnaSequence, possibleStartCodonIndexes3to5, possibleStopCodonIndexes3to5);

        allORFs.addAll(orFs5to3);
        allORFs.addAll(orFs3to5);
        return allORFs;
    }

    /**
     * Find/Create ORFs given sequence and start and stop codon indexes in 5' to 3' direction.
     *      Making sure to prune ORFs that have nested ORFs. IE take longest ORF
     * @param dnaSequence sequence to use
     * @param possibleStartCodonIndexes5to3 list of possible start codon indexes
     * @param possibleStopCodonIndexes5to3  list of possible stop codon indexes
     * @return list of ORFs found
     */
    private static ArrayList<ORF> createORFs5to3(StringBuilder dnaSequence, ArrayList<Integer> possibleStartCodonIndexes5to3, ArrayList<Integer> possibleStopCodonIndexes5to3) {
        ArrayList<ORF> matches = new ArrayList<>();
        if (possibleStartCodonIndexes5to3.isEmpty() || possibleStopCodonIndexes5to3.isEmpty())
            return matches;

        for (Integer startCodonIndex : possibleStartCodonIndexes5to3) {
            ArrayList<Integer> stopCodons = (ArrayList<Integer>) possibleStopCodonIndexes5to3.clone();
            // dont check stop codons that happen before the start codon,
            //      and are in a different reading frame
            stopCodons.removeIf(s ->
                    s <= startCodonIndex ||
                            (startCodonIndex % 3 != s % 3));
            if (!stopCodons.isEmpty()) {
                boolean nestedSequence = false;
                Collections.sort(stopCodons);
                Integer stopCodonIndex = stopCodons.get(0);
                // Check if other ORFs already have been made with this stopindex and skip current ORF if so.
                //      If so, it was a longer ORF then the one currently being processed due to order
                for(ORF orf: matches){
                    if(orf.stopIndex == stopCodonIndex)
                        nestedSequence = true;
                }
                if(!nestedSequence){
                    ORF orf = new ORF(startCodonIndex, stopCodonIndex, dnaSequence.substring(startCodonIndex, stopCodonIndex + 3), 0);
                    ORF.calcORFAcidSeq(orf);
                    matches.add(orf);
                }
            }
        }
        return matches;
    }

    /**
     * Find/Create ORFs given sequence and start and stop codon indexes in 3' to 5' direction.
     *      Making sure to prune ORFs that have nested ORFs. IE take longest ORF
     * @param dnaSequence sequence to use
     * @param possibleStartCodonIndexes3to5 list of possible start codon indexes
     * @param possibleStopCodonIndexes3to5  list of possible stop codon indexes
     * @return list of ORFs found
     */
    private static ArrayList<ORF> createORFs3to5(StringBuilder dnaSequence, ArrayList<Integer> possibleStartCodonIndexes3to5, ArrayList<Integer> possibleStopCodonIndexes3to5){
        ArrayList<ORF> matches = new ArrayList<>();
        /*
            Reverse order of start codons so that we look at largest start codon index first,
                Allows logic to check for nested ORFs to work because of order of matching
         */
        possibleStartCodonIndexes3to5.sort(Collections.reverseOrder());
        if (possibleStartCodonIndexes3to5.isEmpty() || possibleStopCodonIndexes3to5.isEmpty())
            return matches;
        for (Integer startCodonIndex : possibleStartCodonIndexes3to5){
            ArrayList<Integer> stopCodons = (ArrayList<Integer>) possibleStopCodonIndexes3to5.clone();
            // dont check stop codons that happen before the start codon,
            //      and are in a different reading frame
            stopCodons.removeIf(s ->
                    s >= startCodonIndex ||
                            (startCodonIndex % 3 != s % 3));
            if (!stopCodons.isEmpty()) {
                boolean nestedSequence = false;
                Collections.sort(stopCodons);
                Integer stopCodonIndex = stopCodons.get(stopCodons.size()-1);
                // Check if other ORFs already have been made with this stopindex and skip current ORF if so.
                //      If so, it was a longer ORF then the one currently being processed due to order
                for(ORF orf: matches){
                    if(orf.stopIndex == stopCodonIndex){
                        nestedSequence = true;
                    }
                }
                if(!nestedSequence){
                    ORF orf = new ORF(startCodonIndex, stopCodonIndex, dnaSequence.substring(stopCodonIndex, startCodonIndex+3), 1);
                    ORF.calcORFAcidSeq(orf);
                    matches.add(orf);
                }
            }
        }
        return matches;
    }

    /**
     * Find and output all Occuring Indexes of a given sequence
     * @param dnaSequence the sequence to examine
     * @param sequenceToMatch the sequence you want to match
     * @return list of indexes where all sequenceToMatch occurrences were found
     */
    private static ArrayList<Integer> findCodonIndexes(StringBuilder dnaSequence, String sequenceToMatch) {
        ArrayList<Integer> possibleCodonIndexes = new ArrayList<>();
        int res = 0;
        while (res != -1) {
            res = dnaSequence.indexOf(sequenceToMatch, res);
            if (res != -1) {
                possibleCodonIndexes.add(res);
                res++;
            }
        }
        return possibleCodonIndexes;
    }

    private static ArrayList<Integer> findPossibleStopCodonIndexes5to3(StringBuilder dnaSequence) {
        ArrayList<Integer> possibleStopCodonIndexes = new ArrayList<>();
        possibleStopCodonIndexes.addAll(findCodonIndexes(dnaSequence, "TAA"));
        possibleStopCodonIndexes.addAll(findCodonIndexes(dnaSequence, "TAG"));
        possibleStopCodonIndexes.addAll(findCodonIndexes(dnaSequence, "TGA"));
        return possibleStopCodonIndexes;
    }

    private static ArrayList<Integer> findPossibleStopCodonIndexes3to5(StringBuilder dnaSequence) {
        ArrayList<Integer> possibleStopCodonIndexes = new ArrayList<>();
        possibleStopCodonIndexes.addAll(findCodonIndexes(dnaSequence, "TTA"));
        possibleStopCodonIndexes.addAll(findCodonIndexes(dnaSequence, "CTA"));
        possibleStopCodonIndexes.addAll(findCodonIndexes(dnaSequence, "TCA"));
        return possibleStopCodonIndexes;
    }

    /*
        Calculate the Acid Sequence of the given ORF
            Taking into account differnet logic depending on direction.
            if ORF direction is 3' to 5' then reverse the sequence and then use the 3to5 map
     */
    private static void calcORFAcidSeq(ORF orf) {
        StringBuilder str = new StringBuilder();
        StringBuilder dnaSequence = new StringBuilder(orf.sequence);
        if(orf.direction == 1)
            dnaSequence.reverse();
        // For each Codon, look it up in the map corresponding to the direction
        for (int i = 0; i < dnaSequence.length(); i += 3) {
            if(orf.direction == 0){
                str.append(acidMap5to3.get(dnaSequence.substring(i, i + 3)));
                continue;
            }
            str.append(acidMap3to5.get(dnaSequence.substring(i, i + 3)));
        }
        orf.acidSequence = str.toString();
    }

    @Override
    public String toString() {
        if(this.acidSequence == null)
            ORF.calcORFAcidSeq(this);
        String direction = "5' --> 3'";
        if(this.direction == 1)
            return "Start: " + (startIndex+3) + "\tStop: " + (stopIndex+1) + "\nDirection: "
                    + "5' <-- 3'" + "\nSequence: "  + sequence + "\nAcid Sequence: " + this.acidSequence;

        return "Start: " + (startIndex+1) + "\tStop: " + (stopIndex+3) + "\nDirection: "
                + "5' --> 3'" + "\nSequence: "  + sequence + "\nAcid Sequence: " + this.acidSequence;
    }

}
