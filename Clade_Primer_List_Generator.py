from Bio.Blast import NCBIStandalone
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
import MySQLdb
from Bio import Seq
from itertools import product
import csv
import sys

# SYSTEM VARIABLES

MONITOR = ""
HOSTlocal = "localhost"
HOSTremote = "localhost"
DB = "Torrent"
USER = "root"
PASS = <password>
ComputerID = ""

Btag = "GTTTCCAAGTCACCTAGT"  #was originally appending to sequences, but commented out that line of code (Line 183)

###############################################

def SetupTables():
    conn = MySQLdb.connect(host = HOSTremote,
                           user = USER,
                           passwd = PASS,
                           db = DB)
    cursor = conn.cursor()
    
    dropline = ("DROP TABLE IF EXISTS VIR_PrimerList") 
    cursor.execute(dropline)
    
    executeline =("""CREATE TABLE VIR_PrimerList
                       (
                       PrID INTEGER,
                       Posn INTEGER,
                       Tm VARCHAR(5),
                       Sense VARCHAR(25)
                       INDEX (PrID))                     
                       """)
    cursor.execute(executeline)

    cursor.close() 
    conn.commit()
    conn.close()

def extend_ambiguous_dna(seq):
   #return list of all possible sequences given an ambiguous DNA input
   d = Seq.IUPAC.IUPACData.ambiguous_dna_values
   r = []
   for i in product(*[d[j] for j in seq]):
      r.append("".join(i))
   return r 


def BLASTN_v29(outfile, infile, dbfile):
    import os
    import sys
    print "BLASTn (Sequence Read Match)"
    BLASTN_CMDLINE = "/usr/local/ncbi/blast/bin/blastn -out %s -query %s -db %s -task blastn-short -outfmt 5"# -num_alignments -task blastn-short 100000" #took out -task blastn-short as it was too slow and didn't increase alignments
    
    status = os.system(BLASTN_CMDLINE % (outfile, infile, dbfile))


def MakeDB(Inputfile, dbType, DBtitle):
    import os
    import sys
    
    print "Building BLAST database..."
    MAKEDB_CMDLINE = "/usr/local/ncbi/blast/bin/makeblastdb -in %s -out %s -parse_seqids -dbtype %s"
    status = os.system(MAKEDB_CMDLINE % (Inputfile, DBtitle, dbType))


def ReadFastaALN_File(ASF_ALN_infile, PrimerLength, PrimerGap): 
    from Bio import SeqIO

    FASTAoutF = "/users/rwbarrettemac/desktop/ASF_WGS/ASF_PrimerSetF.fasta"
    FASTAoutR = "/users/rwbarrettemac/desktop/ASF_WGS/ASF_PrimerSetR.fasta"

    CSVout = "/users/rwbarrettemac/desktop/ASF_WGS/ASF_PrimerSet.csv"
    CSVoutF = "/users/rwbarrettemac/desktop/ASF_WGS/ASF_PrimerSetF.csv"
    CSVoutR = "/users/rwbarrettemac/desktop/ASF_WGS/ASF_PrimerSetR.csv"


    fasta_fileF = open(FASTAoutF, "w")
    fasta_fileR = open(FASTAoutR, "w")

    #fasta_fileF.write(">" + title + "\n")
    #fasta_fileR.write(">" + title + "\n")

    #for i in range(0, len(sequence), 72):
    #    fasta_file.write(sequence[i:i+72])
    #    fasta_file.write("\n")

    #for i in range(0, len(sequence), 72):
    #    fasta_file.write(sequence[i:i+72])
    #    fasta_file.write("\n")
    #fasta_fileF.close()
    #fasta_fileR.close()

    f = open(CSVout, "w")
    writer = csv.writer(f)
    writer.writerow( ('Primer_Count', 'Primer_ID', 'Primer_Seq', 'Sense', 'Tm', 'ALN_target') )

    ff = open(CSVoutF, "w")
    writerF = csv.writer(ff)
    writerF.writerow( ('Primer_Count', 'Primer_ID', 'Primer_Seq', 'Sense', 'Tm', 'ALN_target') )

    fr = open(CSVoutR, "w")
    writerR = csv.writer(fr)
    writerR.writerow( ('Primer_Count', 'Primer_ID', 'Primer_Seq', 'Sense', 'Tm', 'ALN_target') )

    #PrimerCount, ID, OrderPrimer, CurSense, Tm
    
    ID = 0
    handle = open(ASF_ALN_infile)
    fq_dict2 = SeqIO.parse(handle,"fasta")

    #print "Reading Fastq " + Direction + "..."
    Sense = "FWD"
    PrimerStart = 0
    PrimerCount = 0
    for seqb in fq_dict2:

               
                ##fastq_file = open(fastq_name, "a")
                ContigID = seqb.description
    #            
                #print ContigID
                fasta_sequence = seqb.seq.tostring()

                fasta_LEN = len(fasta_sequence)
                for prSeq in range(0,fasta_LEN):

                    ID = ID + 1

                    PrimerStart = PrimerStart+1
                    PrimerEnd = PrimerStart+PrimerLength
                    PrimerScan = fasta_sequence[PrimerStart:PrimerEnd]


                    aPscan = PrimerScan.count("A")
                    cPscan = PrimerScan.count("C")
                    tPscan = PrimerScan.count("T")
                    gPscan = PrimerScan.count("G")

                    totalPscan = aPscan+cPscan+tPscan+gPscan
                  
                    if Sense == "FWD":
                        #OrderPrimer = tp
                        Sense = "REV"
                        CurSense = "FWD"
                    else:
                        #OrderPrimer = RevPrimer
                        Sense = "FWD"
                        CurSense = "REV"

                    if totalPscan >= PrimerLength-5:  # Ambiguiit
                        testPrimers = extend_ambiguous_dna(PrimerScan)
                        if CurSense == "FWD":
                            PrimerStart = PrimerStart+PrimerLength+PrimerGap# was -12 
                            ID = ID+PrimerLength+PrimerGap #was -12
                        else:
                            PrimerStart = PrimerStart+PrimerLength+PrimerGap #was +1
                            ID = ID+PrimerLength+PrimerGap #was +1
                        for tp in testPrimers:

                            Tm, GCpct = CalcTm(tp)
                            #if Tm >= 51 and Tm <=58:
                            if GCpct >=0 and GCpct <=100:

                                    PrimerCount = PrimerCount+1

                                    ###PrimerStart = PrimerStart+PrimerLength-12
                                    ###ID = ID+PrimerLength-12
                                    #print(PrimerCount, ID, tp, Tm)
                                    from Bio.Seq import Seq

                                    seq = Seq(str(tp))
                                    
                                    RevPrimer = seq.reverse_complement()
                                    #OrderPrimer = Btag+RevPrimer

                                    if Sense == "FWD":
                                        OrderPrimer = str(seq[0:])  
                                        print(OrderPrimer)                           
                                    else:

                                        OrderPrimer = str(RevPrimer[0:])
                                        print(OrderPrimer)

                                    SeqID = str(ID)+CurSense
                                    #print(PrimerCount, ID, OrderPrimer, CurSense, Tm)
                                    ALN_ID = BlastSequencesV29(OrderPrimer, SeqID, PrimerLength)
                                    if ALN_ID != "NA":
                                        # write output to CSV file

                                        #'Primer_Count', 'Primer_ID', 'Primer_Seq', 'Sense', 'Tm', 'ALN_target'
                                       
                                        writer.writerow( (PrimerCount, ID, OrderPrimer, CurSense, Tm, ALN_ID, GCpct ))
                                        print(PrimerCount, ID, OrderPrimer, CurSense, Tm, ALN_ID)
                                        
                                        title = str(PrimerCount)+"_"+str(ID)+"_"+CurSense+str(PrimerCount)
                                        if CurSense == "FWD":
                                            writerF.writerow( (PrimerCount, ID, OrderPrimer, CurSense, Tm, ALN_ID, GCpct ))

                                            fasta_fileF.write(">" + title + "\n")
                                            for i in range(0, len(OrderPrimer), 72):
                                                fasta_fileF.write(OrderPrimer[i:i+72])
                                                fasta_fileF.write("\n")
                                        else:
                                            writerR.writerow( (PrimerCount, ID, OrderPrimer, CurSense, Tm, ALN_ID, GCpct ))

                                            fasta_fileR.write(">" + title + "\n")
                                            for i in range(0, len(OrderPrimer), 72):
                                                fasta_fileR.write(OrderPrimer[i:i+72])
                                                fasta_fileR.write("\n")
                                    

                                    
    fasta_fileF.close()
    fasta_fileR.close()

    f.close()
    ff.close()
    fr.close()




def CalcTm(Oligo):
    import math

# Set initial parameters

    HitOligo = Oligo
    HitOligComp = ""        
    A_countH = 0
    T_countH = 0
    C_countH = 0
    G_countH = 0    
    OligoLen = len(HitOligo)

# Get Complement of BLAST hit Oligo

    for i1 in range(0,OligoLen):
            HitOligAA = HitOligo[i1:i1+1]
            if HitOligAA == "A":
                A_countH = A_countH + 1
                HitOligComp = HitOligComp+"T"
            elif HitOligAA == "T":
                T_countH = T_countH + 1
                HitOligComp = HitOligComp+"A"
            elif HitOligAA == "C":
                C_countH = C_countH + 1
                HitOligComp = HitOligComp+"G"
            elif HitOligAA == "G":
                G_countH = G_countH + 1
                HitOligComp = HitOligComp+"C"

# Calculate Tm for BLAST hit and Subject
# Method by: Sambrook,J., and Russell,D.W. (2001) Molecular Cloning: A Laboratory Manual.
#            Cold Spring Harbor Laboratory Press; Cold Spring Harbor, NY.

    TmAH = 79.8 + 18.5 * math.log10(SaltConc) + (58.4*(G_countH + C_countH)/(A_countH + T_countH + C_countH + G_countH))
    TmBH = (11.8 * ((G_countH + C_countH)/(A_countH + T_countH + C_countH + G_countH))**2)-(820/(A_countH + T_countH + C_countH + G_countH))
    TmH = TmAH + TmBH
    GCratio = (float(G_countH + C_countH)/OligoLen)*100
    #print(GCratio)
    OligoTm = TmH

    #print "Tm (BLAST hit): " + str(round(TmH,2)) + "'C"
    return OligoTm, GCratio

def Save_fasta(filename, title, sequence):
    fasta_file = open(filename, "w")
    fasta_file.write(">" + title + "\n")
    for i in range(0, len(sequence), 72):
        fasta_file.write(sequence[i:i+72])
        fasta_file.write("\n")
    fasta_file.close()


def BlastSequencesV29(SeqIn, SeqID, TargetLEN):

    from Bio.Blast.Applications import NcbiblastnCommandline 

    Save_fasta("/users/rwbarrettemac/bioinformatics/PrimerSeq.fasta", SeqID, SeqIn)

    PreGI = "NA"
    IDcnt = 0
    INDEX = 0
    preIndexCNT = 0
    DBblast = "/users/rwbarrettemac/Bioinformatics/BlastDB/VIRdb"
    #DBblast = "/users/rwbarrettemac/Bioinformatics/BlastDB/NTdb"
    #SaveTitle = "c:\\TorrentFiles\\FastaQue.fasta"
    print "Running BLAST..."
    
    BLASTN_v29("/users/rwbarrettemac/bioinformatics/pythonfolders/TorrentFiles/PrimerCHECK.xml", "/users/rwbarrettemac/bioinformatics/PrimerSeq.fasta", DBblast)
                           
    from Bio.Blast import NCBIXML
    result_handleNEW = open("/users/rwbarrettemac/bioinformatics/pythonfolders/TorrentFiles/PrimerCHECK.xml")
    blast_records = NCBIXML.parse(result_handleNEW)
    for blast_record in blast_records:
            
            
                            for alignment in blast_record.alignments:
                                preIndexCNT = preIndexCNT+1
                                if IDcnt <1:
                                   for hsp in alignment.hsps:
                                       if hsp.identities ==TargetLEN:
                                        IDcnt = IDcnt +1
                                        IDpep = float(hsp.identities)
                                        Lenpep = float(len(hsp.query))
                                        if Lenpep>=TargetLEN: # Had been changing this value, but decided to target it to primer length
                                           prcntid = IDpep/Lenpep
                                              
                                           PreGI = alignment.title
                                           print PreGI, ".......",Lenpep, "nt"
                                           QueSEQ = hsp.query  # FASTA in file
                                           MatchSEQ = hsp.match
                                           SubjSEQ = hsp.sbjct  # Illumina record
                                           MatchCNT = MatchSEQ.count("|")
                                           #print SubjSEQ

                                           
                                           StartPos = hsp.query_start
                                           EndPos = hsp.query_end
    return PreGI        
                                    

        
SaltConc = 0.05
MinTemp = 40
MaxTemp = 90

#MakeDB("/users/rwbarrettemac/Bioinformatics/BlastDB/nt.fasta","nucl", "/users/rwbarrettemac/Bioinformatics/BlastDB/NTdb")


MakeDB("/users/rwbarrettemac/Downloads/ASF_AllCG.fasta", "nucl", "/users/rwbarrettemac/Bioinformatics/BlastDB/VIRdb")


# ReadFastaALN_File( INPUT FILE AS FASTA CONSENSUS SEQUENCE (Single), DESIRED LENGTH OF PRIMER, MINIMAL NT GAP BETWEEN PRIMER SELECTIONS)

# **** See Line 80 -> 85 for output file location ****.

ReadFastaALN_File("/users/rwbarrettemac/desktop/ASF_WGS/ASF_CONS_TEST_FINAL.fasta", 80, 150)
