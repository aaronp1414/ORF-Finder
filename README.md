# 123A_Project (ORF Finder Given DNA Sequence)
  This program was made to implement a fundamental bioinformatics topic involving gene prediction. Specifically this program will take a sequence of DNA as an input, and then output all possible Open Reading Frames (ORFs), including 3' to 5' and 5' to 3' directions.
  
  Open Reading Frames are the 'code' used for creating proteins. Knowing the possible Open Reading Frames is crucial in finding spots in DNA that result in certain protens and functions in our bodies. 
  
  ## Input
  Given a File containing a DNA sequence, preferably FASTA format, this program will parse through the file and output the open reading frames that could result from the given DNA sequence
  
  ## Output
  Output will be written to console. If ORFs are nested, the longest ORF will be kept and nested ORFs will be ignored. 
  
  https://www.ncbi.nlm.nih.gov/orffinder/ was used to check this program for accuracy and found no differences in found ORFs
  
  ## Sample Output
  Sample Output:
  
  ![Example Output](https://i.gyazo.com/9bd42ee623e835e1c72578235b4f61a7.png)
