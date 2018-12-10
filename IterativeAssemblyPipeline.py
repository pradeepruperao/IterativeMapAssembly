#!/usr/bin/env python
# Author: Pradeep Ruperao 
# Usage: python IterativeAssembly.py -l pyfqlist.txt -c pyconv.txt -r /main/projects/cassava/GenomeComparison/graph/Ref/NC_003070.fa

import os,sys
from optparse import OptionParser
from shutil import copyfile
import subprocess

fastq={}
conv={}
iter=1



def getdata(fastqList,conventionFileName,reference):
  fp = open(conventionFileName, 'r')
  for line in fp:
    cn=[]
    cn.append(line.rstrip('\n').split('\t'))
    conv[cn[0][0]]=cn[0][1]
  fq=open(fastqList,'r')
  for fstq in fq:
    cn=[]
    cn.append(fstq.rstrip('\n').split('\t'))
    fastq[cn[0][0]]=cn[0][1]

def runIterAsbly(reference):
  novelseq=""
  pwd=os.getcwd()
  for con in conv:
    outfile = open("cmd.sh","w")
    accdir=pwd+"/"+conv[con]
    refdir=accdir+"/01_Ref"
    if not os.path.exists(accdir):
      os.mkdir(accdir)
    if not os.path.exists(refdir):
      os.mkdir(refdir)
    refname=os.path.basename(reference)
    refname=conv[con]+".fa"
    if novelseq != "":
      cmd="cd "+refdir+"\n cat "+ reference+" "+novelseq+">"+refname
      #os.system(cmd)
      outfile.write("\n"+cmd) 
    else:
      copyfile(reference, refdir+"/"+refname)
    cmd="  bowtie2-build " + refdir + "/" + refname + " " + refdir + "/" + refname
    #os.system(cmd); # Index reference file
    outfile.write("\n"+cmd)
    
    #Trim fastq data
    trimdir=accdir+"/02_Trimming"
    if not os.path.exists(trimdir):
      os.mkdir(trimdir)
    for file in fastq:
      if con in file:
        cmd="\ncd "+trimdir+"\n"
        outfile.write("\n"+cmd)
        R1=os.path.basename(file)
        R2=os.path.basename(fastq[file])
        #print(file,fastq[file])
        trimcmd="java -jar /scratch/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 2 -phred33 " + file + " " + fastq[file] +" "+ trimdir+"/"+ R1 +".fwd.pair.gz"+" "+ trimdir +"/"+ R1 +".fwd.single.gz"+" "+ trimdir+"/"+ R2 +".rev.pair.gz"+" "+ trimdir+"/"+ R2 +".rev.single.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:10:25 MINLEN:35 >>Trim.log 2>&1 ";
        #os.system(trimcmd)
        outfile.write("\n"+trimcmd)
        
    #Map reads
    mapdir=accdir+"/03_Mapping"
    if not os.path.exists(mapdir):
      os.mkdir(mapdir)
    for file in fastq:
      if con in file:
        R1=os.path.basename(file)
        R2=os.path.basename(fastq[file])
        cmd="\n cd "+mapdir+"\n"
        outfile.write("\n"+cmd)
        out=R1+".sort.bam";
        mapcmd="bowtie2 -1 "+trimdir+"/"+R1+".fwd.pair.gz -2 "+trimdir+"/"+R2+".rev.pair.gz -x "+refdir+"/"+refname+" -p 10 -q --end-to-end --sensitive -I 0 -X 1000 | samtools view -Sb - | samtools sort - -o "+mapdir+"/"+out+" &>Pipeline.log ";
        #os.system(mapcmd)
        outfile.write("\n"+mapcmd)
        idxcmd="samtools index "+mapdir+"/"+out
        #os.system(idxcmd)
        outfile.write("\n"+idxcmd)
        idslistCmd="samtools view -f 12 -F 256 "+mapdir+"/"+out+" |cut -f 1|sort|uniq >"+mapdir+"/temp.txt";
        #os.system(idslistCmd)
        outfile.write("\n"+idslistCmd)
        rawfq1cmd="perl ~/scripts/retrieve_fastq_fast.pl -in "+file+" -out "+mapdir+"/R1.fastq -l "+mapdir+"/temp.txt";
        rawfq2cmd="perl ~/scripts/retrieve_fastq_fast.pl -in "+file+" -out "+mapdir+"/R2.fastq -l "+mapdir+"/temp.txt";
        #os.system(rawfq1cmd)
        #os.system(rawfq2cmd)
        outfile.write("\n"+rawfq1cmd)
        outfile.write("\n"+rawfq2cmd)
        
    # Assembly
    ablydir=accdir+"/04_Assembly";
    if not os.path.exists(ablydir):
      os.mkdir(ablydir)    
    R1=os.path.basename(file)
    asblycmd="spades.py -t 20  -m 800 -1 "+mapdir+"/R1.fastq -2 "+mapdir+"/R2.fastq -o "+ablydir+"/"+R1+"Spadesout"
    cmd="\n cd "+ablydir+"\n"
    outfile.write("\n"+cmd)
    outfile.write("\n"+asblycmd)
    novelseq=ablydir+"/"+R1+"Spadesout/scaffolds.fasta";
    outfile.close()
    #break
    os.system("bash cmd.sh")




if __name__ == '__main__':
  parser = OptionParser("Usage: %prog -f  fastq files [Option]")

  parser.add_option("-l", "--fql", type="string", dest="fastqList",
                    help="Give a list of raw R1 read files ")

  parser.add_option("-c", "--cnv", type="string", dest="conventionFileName",
                   help="Give a tab-delimited convention file name. EG: A1	AccName")  

  parser.add_option("-r", "--reference", type="string", dest="reference", 
                    help="Give reference to start assembly")

  (opts, args) = parser.parse_args()

  if (opts.fastqList is None):
    print("Give a list of raw R1 reads files")
    parser.print_help()
  elif (opts.conventionFileName is None):
    print ("Give a tab-delimited convention file name. EG: A1	AccName")
    parser.print_help()
  elif(opts.reference is None):
    print("Give reference sequence file")
  else:
    getdata(opts.fastqList,opts.conventionFileName,opts.reference)
    runIterAsbly(opts.reference)



