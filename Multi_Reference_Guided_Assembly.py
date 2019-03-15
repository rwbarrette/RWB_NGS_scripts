from Bio.Blast import NCBIStandalone
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
import MySQLdb
import os
import string


        
def RunVelvet(InfileFWD, InfileREV):
    import os
    import sys
    print "Running Velveth"

    VELVETh_CMDLINE = "/users/rwbarrettemac/bioinformatics/velvet/velveth /users/rwbarrettemac/bioinformatics/velvet/data/assem/FWDreads 21 -fastq -short %s"
    status = os.system(VELVETh_CMDLINE % (InfileFWD))

    VELVETg_CMDLINE = "/users/rwbarrettemac/bioinformatics/velvet/velvetg /users/rwbarrettemac/bioinformatics/velvet/data/assem/FWDreads -unused_reads yes"
    status = os.system(VELVETg_CMDLINE)

    VELVETh_CMDLINE = "/users/rwbarrettemac/bioinformatics/velvet/velveth /users/rwbarrettemac/bioinformatics/velvet/data/assem/REVreads 21 -fastq -short %s"
    status = os.system(VELVETh_CMDLINE % (InfileREV))

    VELVETg_CMDLINE = "/users/rwbarrettemac/bioinformatics/velvet/velvetg /users/rwbarrettemac/bioinformatics/velvet/data/assem/REVreads -unused_reads yes"
    status = os.system(VELVETg_CMDLINE)

# SYSTEM VARIABLES

MONITOR = "roger.w.barrette@gmail.com"
HOSTlocal = "localhost"
HOSTremote = "localhost"
DB = "Torrent"
USER = <user>
PASS = <password>
ComputerID = "NODE1"


#RunVelvet("/users/rwbarrettemac/bioinformatics/BWAfiles/FMDsat3v3529_S5_L001_R1_001.fastq","/users/rwbarrettemac/bioinformatics/BWAfiles/FMDsat3v3529_S5_L001_R2_001.fastq")


def RunBWA(REFInfile, FastqInFileFWD, FastqInFileREV):
    import os
    import sys
    print "Running BWA Indexing and SAI file build..."
    

    BWA_index = "/users/rwbarrettemac/bioinformatics/bwa/bwa index %s" % (REFInfile)
    status = os.system(BWA_index)

    BWA_saiF = "/users/rwbarrettemac/bioinformatics/bwa/bwa mem %s %s %s >  /users/rwbarrettemac/bioinformatics/BWAfiles/BWAsam.sam" % (REFInfile, FastqInFileFWD, FastqInFileREV)
    status = os.system(BWA_saiF)

def RunBWA_SE(REFInfile, FastqInFileFWD, FastqInFileREV):
    import os
    import sys
    print "Running BWA Indexing and SAI file build..."
    
    
    BWA_index = "/users/rwbarrettemac/bioinformatics/bwa/bwa index %s" % (REFInfile)
    status = os.system(BWA_index)
    print BWA_index
    
    BWA_saiF = "/users/rwbarrettemac/bioinformatics/bwa/bwa mem %s %s >  /users/rwbarrettemac/bioinformatics/BWAfiles/BWAsam.sam" % (REFInfile, FastqInFileFWD)
    print BWA_saiF
    
    status = os.system(BWA_saiF)

def RunSAM_SortedBam(REFInfile):
    import os
    import sys
    print "Running SAMTOOLS Index Build (faidx)..."
    
    RunSamIDX = "/users/rwbarrettemac/bioinformatics/samtools/samtools faidx %s" % (REFInfile)
    #print RunSamIDX, "1"
    status = os.system(RunSamIDX)
  
    RunSam_Import = "/users/rwbarrettemac/bioinformatics/samtools/samtools import %s.fai /users/rwbarrettemac/bioinformatics/BWAfiles/BWAsam.sam /users/rwbarrettemac/bioinformatics/BWAfiles/BWAbam.bam" %(REFInfile)
    #print RunSam_Import, "2"
    status = os.system(RunSam_Import)
    
    RunSamSORT = "/users/rwbarrettemac/bioinformatics/samtools/samtools sort /users/rwbarrettemac/bioinformatics/BWAfiles/BWAbam.bam /users/rwbarrettemac/bioinformatics/BWAfiles/BWAbam.sorted"
    #print RunSamSORT, "3"
    status = os.system(RunSamSORT)
    

    RunSamIndexBam = "/users/rwbarrettemac/bioinformatics/samtools/samtools index /users/rwbarrettemac/bioinformatics/BWAfiles/BWAbam.sorted.bam"
    #print RunSamIndexBam, "4"
    status = os.system(RunSamIndexBam)

def RunBAM_STATS(OutfileSTAT):
    import os
    import sys
    print "Running SAMTOOLS Flagstat for read statistics..."
    
    RunSTATS_CMDLINE = "/users/rwbarrettemac/bioinformatics/samtools/samtools flagstat /users/rwbarrettemac/bioinformatics/BWAfiles/BWAbam.bam > %s" %(OutfileSTAT)
    
    status = os.system(RunSTATS_CMDLINE)



def RunSAM_Consensus(REFInfile, Outfile):
    import os
    import sys
    print "Running SAMTOOLS Consensus Build (mpileup)..."
    
    RunSam_Cons_CMDLINE = "/users/rwbarrettemac/bioinformatics/samtools/samtools mpileup -uf %s /users/rwbarrettemac/bioinformatics/BWAfiles/BWAbam.sorted.bam | /users/rwbarrettemac/bioinformatics/samtools/bcftools/bcftools view -cg - | /users/rwbarrettemac/bioinformatics/samtools/bcftools/vcfutils.pl vcf2fq > %s" % (REFInfile, Outfile)
    print RunSam_Cons_CMDLINE
    
    status = os.system(RunSam_Cons_CMDLINE)

def Convert_CONSENS_to_FASTA(Outfile, NewOutfile, OutfileSTAT):
    StatOUTfile = open(OutfileSTAT, "a")
    ModOUTFILE = open(NewOutfile, "w")
    ModINFILE = open(Outfile, "r")
    mLineTotal = 0
    mLineNonN = 0
    inlineCNT = 0
    for mLine in ModINFILE:
        inlineCNT=inlineCNT+1
        if "+" in mLine:
            "DONE"
            break
        else:
            if "@" in mLine:
                StartString = ">"+mLine
                ModOUTFILE.write(StartString)
            
            else:
                #print(mLine)
                mLineLEN = len(mLine)-1
                mLineNcnt = mLine.count("n")
                #print(mLineLEN, mLineNcnt, "**********************###############")
                ModOUTFILE.write(mLine)
                mLineTotal = mLineTotal+(mLineLEN)
                mLineNonN = mLineNonN+(mLineLEN-mLineNcnt)
                #print(mLineTotal, mLineNonN)
                prcntid = (float(mLineNonN/float(mLineTotal)))*100
    
    StatOUTfile.write("\n")
    StatOUTfile.write("### Total NT: " + str(mLineTotal))
    StatOUTfile.write("\n")
    StatOUTfile.write("### ALN NT: " + str(mLineNonN))
    StatOUTfile.write("\n")
    StatOUTfile.write("### Percent ID(%): "+ str(prcntid))
    StatOUTfile.write("\n")


    ModOUTFILE.close()
    ModINFILE.close()
    StatOUTfile.close()

    return mLineTotal, mLineNonN     

def ReadFastQinFile(infile_FASTq, outfile_FASTa, New): #Sets up FASTA file from FASTQ
    from Bio import SeqIO
    
    conn = MySQLdb.connect(host = HOSTremote,
                           user = USER,
                           passwd = PASS,
                           db = DB)
    cursor = conn.cursor()
                           
    ID = 0
    seqcount= 0
    seqcountMAX = 100000000
    handleFWD = open(infile_FASTq)
                           
    fq_dictFWD = SeqIO.parse(handleFWD,"fastq")
    if New == "TRUE":

        fastq_file = open(outfile_FASTa, "w")
    else:
        fastq_file = open(outfile_FASTa, "a")

    print "Reading Fastq..."
                           
                           
    for seqFWD in fq_dictFWD:
                        if seqcount <=seqcountMAX:   
                               seqcount = seqcount + 1
                               
                               QscoreMin = min(seqFWD.letter_annotations["phred_quality"])
                               #print QscoreMin
                               ID = ID + 1
                               
                               fastq_sequence = seqFWD.seq.tostring()
                               SeqTitle = seqFWD.description
                               
                               if len(fastq_sequence) >= 0:
                                   
                                   fastq_IDstring = ">"+ SeqTitle + "\n"
                                   fastq_file.write(fastq_IDstring)
                                   
                                   #ArguementsDB = " '"+ str(ID) + "','" +str(fastq_IDstring) + "','" + fastq_sequence + "','" + str(QscoreMin) + "','" + str("FWD") + "' "
                                   #EnterLine = "fastQReads(SeqGI, SeqID, Sequence, Q_score, ReadFR) VALUES (" + ArguementsDB + ")"
                                   #ActLine = "INSERT INTO " + EnterLine
                                   
                                   #cursor.execute(ActLine)
                                   
                                   for i in range(0, len(fastq_sequence), 72):
                                       fastq_file.write(fastq_sequence[i:i+72])
                                       fastq_file.write("\n")
    fastq_file.close()
    handleFWD.close()
    
    cursor.close()
    conn.commit()
    conn.close()

def RunBowtie(RefInfile, FastqInFileFWD, FastqInFileREV):

    import os
    import sys
    print "Running Bowtie Indexing and SAI file build..."
    
    Bow_index = "/users/rwbarrettemac/bioinformatics/bowtie/bowtie-build %s /users/rwbarrettemac/bioinformatics/bowtie/bowtieREF" % (RefInfile)
    status = os.system(Bow_index)

    Bow_process = "/users/rwbarrettemac/bioinformatics/bowtie/bowtie -S -p 2 /users/rwbarrettemac/bioinformatics/bowtie/bowtieREF -1 %s -2 %s /users/rwbarrettemac/bioinformatics/BWAfiles/BWAsam.sam" % (FastqInFileFWD, FastqInFileREV)
    status = os.system(Bow_process)

def RunBowtie2(RefInfile, FastqInFileFWD, FastqInFileREV):

    import os
    import sys
    print "Running Bowtie2 Indexing and SAI file build..."

    Bow2_index = "/users/rwbarrettemac/bioinformatics/bowtie2/bowtie2-build %s /users/rwbarrettemac/bioinformatics/bowtie2/bowtieREF" % (RefInfile)
    status = os.system(Bow2_index)
    
    Bow2_process = "/users/rwbarrettemac/bioinformatics/bowtie2/bowtie2 -x /users/rwbarrettemac/bioinformatics/bowtie2/bowtieREF -1 %s -2 %s -S /users/rwbarrettemac/bioinformatics/BWAfiles/BWAsam.sam --un-conc /users/rwbarrettemac/bioinformatics/bowtie2/unmapped.fastq --al-conc /users/rwbarrettemac/bioinformatics/bowtie2/mapped.fastq" % (FastqInFileFWD, FastqInFileREV)
    status = os.system(Bow2_process)

def RunGsnap(REFInfile, FastqInFileFWD, FastqInFileREV):
    import os
    import sys
    print "Running GSnap Indexing and SAI file build..."
    

    GSnap_index = "gmap_build -d ViralRefGMp %s" % (REFInfile)
    print GSnap_index
    status = os.system(GSnap_index)

    GSnap_saiF = "gsnap -d ViralRefGMp %s %s -A sam > /users/rwbarrettemac/bioinformatics/BWAfiles/BWAsam.sam" % (FastqInFileFWD, FastqInFileREV)
    #GSnap_saiF = "gsnap -d ViralRefGMp %s -A sam > /users/rwbarrettemac/bioinformatics/BWAfiles/BWAsam.sam" % (FastqInFileFWD)
    status = os.system(GSnap_saiF)

def BLASTN_v29(outfile, infile, dbfile):
    import os
    import sys
    print "BLASTn (Sequence Read Match)"
    BLASTN_CMDLINE = "/usr/local/ncbi/blast/bin/blastn -out %s -query %s -db %s -outfmt 5 -num_alignments 1000000"# -num_alignments -task blastn-short 100000" #took out -task blastn-short as it was too slow and didn't increase alignments
    
    status = os.system(BLASTN_CMDLINE % (outfile, infile, dbfile))


def BlastSequencesV29(TemplateInFile, BlastType, OutFileStat):

    from Bio.Blast.Applications import NcbiblastnCommandline
    
    
    
    BlastNThandle = open(TemplateInFile)
    BlastNTread = SeqIO.parse(BlastNThandle,"fasta")
    for reads in BlastNTread:
        TestSeq = reads.seq.tostring()
    
    QueLen = len(TestSeq)
    print "Template Length: ", QueLen
    
    IDcnt = 0
    INDEX = 0
    preIndexCNT = 0
    DBblast = "/users/rwbarrettemac/bioinformatics/pythonfolders/TorrentFiles/TorrentDB/RawSequence_DIRECT"
    #SaveTitle = "c:\\TorrentFiles\\FastaQue.fasta"
    print "Running BLAST..."
    
    BLASTN_v29("/users/rwbarrettemac/bioinformatics/pythonfolders/TorrentFiles/VIRCurrentTorrent_DIR.xml", TemplateInFile, DBblast)
    
    fastq_file = open(OutFileStat, "a")

    print "BLAST complete."
    
    conn = MySQLdb.connect(host = HOSTlocal,
                           user = USER,
                           passwd = PASS,
                           db = DB)
    cursor = conn.cursor()
                           
    from Bio.Blast import NCBIXML
    result_handleNEW = open("/users/rwbarrettemac/bioinformatics/pythonfolders/TorrentFiles/VIRCurrentTorrent_DIR.xml")
    blast_records = NCBIXML.parse(result_handleNEW)
    for blast_record in blast_records:
            
            
                               for alignment in blast_record.alignments:
                                   preIndexCNT = preIndexCNT+1
                #if IDcnt <1:
                                   for hsp in alignment.hsps:
                                       if hsp.identities >=0:
                                        
                                        IDpep = float(hsp.identities)
                                        Lenpep = float(len(hsp.query))
                                        if Lenpep>=0:
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
                        
                                           GIstart = 3
                                           PreGIend = PreGI[GIstart:15]
                                           GIend = PreGIend.find("|")+ 3
                                           GIid = PreGI[GIstart:GIend]
                                           INDEX = INDEX + 1

                                           fastq_file.write("\n")
                                           fastq_file.write(PreGI)
                                           fastq_file.write("### Percent ID(%): "+ str(prcntid))
                                           fastq_file.write("\n")
        
    fastq_file.close()
    print "preIndex Count",preIndexCNT
    print "index",INDEX
       


def RunGsnapSE(REFInfile, FastqInFileFWD):
    import os
    import sys
    print "Running GSnap Indexing and SAI file build..."
    

    GSnap_index = "gmap_build -d ViralRefGMp %s" % (REFInfile)
    print GSnap_index
    status = os.system(GSnap_index)

    #GSnap_saiF = "gsnap -d ViralRefGMp %s %s -A sam > /users/rwbarrettemac/bioinformatics/BWAfiles/BWAsam.sam" % (FastqInFileFWD, FastqInFileREV)
    GSnap_saiF = "gsnap -d ViralRefGMp %s -A sam > /users/rwbarrettemac/bioinformatics/BWAfiles/BWAsam.sam" % (FastqInFileFWD)
    status = os.system(GSnap_saiF)

def BWA_SAM(REFInfile,FastqInFileFWD,FastqInFileREV,OutfileFASTa, OutfileSTAT):
    
    FastqInFileFWDassem = FastqInFileFWD #"/users/rwbarrettemac/bioinformatics/velvet/data/assem/FWDreads/contigs.fa"
    FastqInFileREVassem = FastqInFileREV #"/users/rwbarrettemac/bioinformatics/velvet/data/assem/REVreads/contigs.fa"
    OutfileFASTq = "/users/rwbarrettemac/bioinformatics/BWAfiles/FastqQUEOUT.fasta"
    #OutfileSTAT = OutfileFASTa+"_STATS_OUT.txt"
    
    #RunVelvet(FastqInFileFWD,FastqInFileREV) # Editing this line to pass thru less variables !!!
    RunBWA_SE(REFInfile, FastqInFileFWDassem, FastqInFileREVassem) #Changed RunBWA_SE to RunBWA 1/15/2015
    #RunBowtie(REFInfile, FastqInFileFWDassem, FastqInFileREVassem)
    #RunBWA(REFInfile, FastqInFileFWDassem, FastqInFileREVassem)
    #RunGsnap(REFInfile, FastqInFileFWDassem, FastqInFileREVassem)
    #RunGsnapSE(REFInfile, FastqInFileFWDassem)
    #RunBowtie2(REFInfile, FastqInFileFWDassem, FastqInFileREVassem)

    RunSAM_SortedBam(REFInfile)
    #RunBAM_STATS(OutfileSTAT) #Added 2/15/2019 RWB

    RunSAM_Consensus(REFInfile, OutfileFASTa)
    RunBAM_STATS(OutfileSTAT)
    mLineTotal, mLineNonN = Convert_CONSENS_to_FASTA(OutfileFASTa, "/users/rwbarrettemac/bioinformatics/BWAfiles/FastqQUE_TEST.fasta", OutfileSTAT)

    #ReadFastQinFile(FastqInFileFWDassem, OutfileFASTq, "TRUE")
    #ReadFastQinFile(FastqInFileREVassem, OutfileFASTq, "FALSE")
    ###FIX_  MakeDB(REFInfile,"nucl","/users/rwbarrettemac/bioinformatics/pythonfolders/TorrentFiles/TorrentDB/RawSequence_DIRECT")
    ###FIX_  BlastSequencesV29("/users/rwbarrettemac/bioinformatics/BWAfiles/FastqQUE_TEST.fasta", "blastn", OutfileSTAT)
    print "done"


def MakeDB(Inputfile, dbType, DBtitle):
    import os
    import sys
    
    print "Building BLAST database..."
    MAKEDB_CMDLINE = "/usr/local/ncbi/blast/bin/makeblastdb -in %s -out %s -parse_seqids -dbtype %s"
    status = os.system(MAKEDB_CMDLINE % (Inputfile, DBtitle, dbType))


def File_Iterate(FilePath, RefGenome):
    path = FilePath
    listing = os.listdir(path)

    for infile in listing:
        try:

          if ".fastq.gz" in infile:
            continue #print("negatory...")
          else:
            if "R2_001" in infile:
              continue#print("negatory")
            else:
              if "._" in infile:
                continue#print("negatory")
              else:
                print("current file is: " + path + infile)
                #print("POSITORY...")
                RunFileR2 = infile.replace('_R1_001.', '_R2_001.')
                RunFileR2x = path+RunFileR2
                RunFileR1x = path+infile
                FileRootID = infile.replace('L001_R1_001.fastq','')
                OutStat = path+FileRootID+"OUT_STAT.txt"
                OutFASTA = path+FileRootID+"OUT_CONS.fasta"
                print(RunFileR1x, RunFileR2x, OutStat, OutFASTA)

                BWA_SAM(RefGenome,RunFileR1x,RunFileR2x,OutFASTA, OutStat)

          #SystemLine = ("/usr/local/bin/pythonw2 /Volumes/RWB_VIRAL/VirImager_Scripts/Build_NN_Dataset_v5_2CH.py -f %s -i %s -a %s -n %s -c %s -p %s")% (infile, FASTAin, CSVout, CSVoutNT, CSVoutCONS, PDFout)
          #print(SystemLine)
          #os.system(SystemLine) 
        except:
          print("COULD NOT PRODUCE AN ASSEMBLY!!!!")
          continue

File_Iterate("/Volumes/RWB_DATA/MinION_3142019/", "/Volumes/RWB_DATA/ASF/VODFeb2019/ASFV_MiSeq02152019/ASF_Georgia.fasta")


#DONT TOUCH ANYTHING ABOVE THIS LINE _______except RWB...And he shouldn't either ... _________________________________________________

# USAGE:: BWA_SAM( Reference sequence file (FASTA) , Forward Reads from Illumina (FASTQ, Reverse Reads from Illumina (FASTQ),  Output FASTA file )

#BWA_SAM("/users/rwbarrettemac/desktop/BWA_Templates_10_19_2017/S11_FMDV_O_Israel.fasta","/users/rwbarrettemac/desktop/FMD_VB_500/11-11_S11_L001_R1_001.fastq", "/users/rwbarrettemac/desktop/FMD_VB_500/11-11_S11_L001_R2_001.fastq", "/users/rwbarrettemac/desktop/FMD_VB_500/S11_FMD.fasta")

#BWA_SAM("/users/rwbarrettemac/MiSeqRuns/VB_ProblemChild/FMD_C5.fasta","/users/rwbarrettemac/MiSeqRuns/VB_ProblemChild/Sing_C3Chall_S1_L001_R1_001.fastq", "/users/rwbarrettemac/MiSeqRuns/VB_ProblemChild/Sing_C3Chall_S1_L001_R2_001.fastq", "/users/rwbarrettemac/MiSeqRuns/VB_ProblemChild/Sing_Chall_C5_BWA.fasta")

#BWA_SAM("/Volumes/RWB_DATA/ASF/VODFeb2019/ASFV_MiSeq02152019/ASF_Georgia.fasta","/Volumes/RWB_DATA/ASF/VODFeb2019/ASFV_MiSeq02152019/S24Ranctrol_S24_L001_R1_001.fastq", "/Volumes/RWB_DATA/ASF/VODFeb2019/ASFV_MiSeq02152019/S24Ranctrol_S24_L001_R2_001.fastq", "/Volumes/RWB_DATA/ASF/VODFeb2019/ASFV_MiSeq02152019/S24DirectCTRL_S24_BWA.fasta")

#BWA_SAM("/users/rwbarrettemac/bioinformatics/BWAfiles/FMDSAT1ref.fasta","/users/rwbarrettemac/bioinformatics/BWAfiles/FMDsat3v3529_S5_L001_R1_001.fastq","/users/rwbarrettemac/bioinformatics/BWAfiles/FMDsat3v3529_S5_L001_R1_001.fastq", "/users/rwbarrettemac/bioinformatics/BWAfiles/FMDsat3OUT.fastq","/users/rwbarrettemac/bioinformatics/BWAfiles/FMDSat3OUT.fasta")
