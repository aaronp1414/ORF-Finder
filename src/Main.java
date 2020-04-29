import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

public class Main {

    public static void main(String[] args) {
        if(args.length < 1){
            System.out.println("Please give name of fasta file as program argument. IE: test.fasta");
        }
        StringBuilder dnaSequence = getDNAFromFile(args[0]);
        ArrayList<ORF> orfs = ORF.findORFs(dnaSequence);
        if(orfs == null || orfs.isEmpty()){
            System.out.println("No ORFs Found");
            System.exit(0);
        }
        for(ORF orf: orfs){
            System.out.println(orf + "\n");
        }
}

    private static StringBuilder getDNAFromFile(String file) {
        StringBuilder sequence = new StringBuilder();
        try {
            Scanner sc = new Scanner(new File(file));
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


    private static Scanner getScannerFromFile(String file) {
        Scanner sc = null;
            try {
                sc = new Scanner(new File(file));
            } catch (FileNotFoundException e) {
                System.out.println("File Not Found. Please try again");
            }
            return sc;
        }

    }
