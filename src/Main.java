import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

public class Main {

    public static void main(String[] args) {
        boolean validInput = false;
        System.out.println("Please enter path/name of fasta file containing one DNA sequence: ");
        Scanner sc = new Scanner(System.in);
        File file = null;
        while(!validInput){
            file = new File(sc.nextLine());
            validInput = file.exists();
            if(!validInput)
                System.out.println("File " + file.getName() + " was not found. Please enter a new file ");
        }
        System.out.println("Please enter minimum ORF length in nucleotides: ");
        int minSeqLength = Integer.parseInt(sc.nextLine());
        sc.close();
        StringBuilder dnaSequence = getDNAFromFile(file);
        ArrayList<ORF> orfs = ORF.findORFs(dnaSequence, minSeqLength);
        if(orfs == null || orfs.isEmpty()){
            System.out.println("No ORFs Found");
            System.exit(0);
        }
        for(ORF orf: orfs){
            System.out.println(orf + "\n");
        }
}

    private static StringBuilder getDNAFromFile(File file) {
        StringBuilder sequence = new StringBuilder();
        try {
            Scanner sc = new Scanner(file);
            String title = sc.nextLine().substring(1);
            System.out.println("Showing ORFs found in " + title);
            while(sc.hasNextLine()){
                sequence.append(sc.nextLine());
            }
        } catch (IOException e) {
            System.out.println("File not found, Please try again");
            System.exit(1);
        }
        return sequence;
    }
}
