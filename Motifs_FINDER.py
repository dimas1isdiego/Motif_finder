###DIEGO CALZADA FRAILE. 31/1/2017

#Go to directory in the terminal and open as: python3 Motif_finder.py InputFileinFASTA.txt OutputFile.txt

#Function: Opening the fasta sequence
def opensequence(MyFile):
  File = open(MyFile,"r")  #Object to connect with the file
  Sequences = []
  IDs = []
  for line in File:
    if not(">" in line):
        sequence = ""
        sequence = sequence + line.upper().strip()
        Sequences.append(sequence)
    else:
        name = line.strip()
        IDs.append(name)
  fileinfo = dict(zip(IDs, Sequences))
  File.close()
  return(fileinfo)

#Function: search all
def SearchAll(regex, seq, pos, emptylist):
  Res = regex.search(seq, pos)
  if Res == None:
    return()
  else:
    emptylist.append(Res.group())
    SearchAll(regex, seq, Res.start()+1, emptylist)
  return(emptylist)
    

#Function: looking for ExoMotifs
def ExoMotif(dic):
    dicmot = {}
    diccount = {}
    RE1 = r"[GUC]G[ACG][GC]"
    RE2 = r"[CUG]CC[UGA]"
    regex1 = re.compile(RE1)
    regex2 = re.compile(RE2)
    for n in dic:
        Exomotif1 = SearchAll(regex1, dic[n],0,[])
        Exomotif2 = SearchAll(regex2, dic[n],0,[])
        if len(Exomotif1) < 1:
          dicmot[n] = Exomotif2
        elif len(Exomotif2) < 1:
          dicmot[n] = Exomotif1
        else:
          dicmot[n] = Exomotif1 + Exomotif2
        diccount[n] = len(dicmot[n])
    return(dicmot,diccount)
                                                
#Function: looking for CLMotifs
def CLMotif(dicti):
    dictmot = {}
    dictcount = {}
    RE1 = r"[AU][CA]CAU[UCG]"
    RE2 = r"A[GA]G[UG][AG]G[UA]A"
    RE3 = r"[CU]U[GU][CA]A[CU][UAG][GAU]"
    regex1 = re.compile(RE1)
    regex2 = re.compile(RE2)
    regex3 = re.compile(RE3)
                                                
    for n in dicti:
        CLmotif1 = SearchAll(regex1, dicti[n],0,[])
        CLmotif2 = SearchAll(regex2, dicti[n],0,[])
        CLmotif3 = SearchAll(regex3, dicti[n],0,[])

        if len(CLmotif1) < 1:
          if len(CLmotif2) < 1:
            if len(CLmotif3) < 1:
              dictmot[n] = []
            else:
              dictmot[n] = CLmotif3
          else:
            if len(CLmotif3) >= 1:
              dictmot[n] = CLmotif2 + CLmotif3
            else:
              dictmot[n] = CLmotif2
              
        elif len(CLmotif2) < 1:
          if len(CLmotif3) < 1:
            dictmot[n] = CLmotif1
          else:
            dictmot[n] = CLmotif1 + CLmotif3
            
        else:
          if len(CLmotif3) < 1:
            dictmot[n] = CLmotif1 + CLmotif2
          else:
            dictmot[n] = CLmotif1 + CLmotif2 + CLmotif3
          
        dictcount[n] = len(dictmot[n])
    return(dictmot,dictcount)

#Function
def hEXOmotifs(dit):
  ditmot={}
  ditcount={}
  RE=r"[GAU][GUA][GAU][CAG][UA][GC]"
  regex = re.compile(RE)

  for n in dit:
    hEXOmotif = SearchAll(regex, dit[n],0,[])
    if len(hEXOmotif) < 1:
      ditmot[n] = []
    else:
      ditmot[n] = hEXOmotif

    ditcount[n] = len(ditmot[n])
  return(ditmot, ditcount)

    
#Function: Saving into a file
def savinginfo(A1, A2, B1, B2, C1, C2):
    MyFile = open(sys.argv[2], "w")
    MyFile.write("EXO-MOTIF AND CL-MOTIF FINDER \n\n")
    MyFile.write("EXOMOTIFS")
    MyFile.write("\n")
    for n in A1:
        MyFile.write(n)
        MyFile.write("\t\t")
        MyFile.write("|")
        MyFile.write("\t")
        for a in A1[n]:
            MyFile.write(a)
            MyFile.write("\t")
        MyFile.write("\n")

    MyFile.write("\n\n")        
    MyFile.write("Counts of ExoMotifs")        
    MyFile.write("\n")
    for n in A2:
        MyFile.write(n)
        MyFile.write("\t\t")
        MyFile.write("|")
        MyFile.write("\t")
        MyFile.write(str(A2[n]))
        MyFile.write("\n")
    MyFile.write("\n\n\n\n")         
    MyFile.write("CLMOTIFS")
    MyFile.write("\n")
    for n in B1:
        MyFile.write(n)
        MyFile.write("\t\t")
        MyFile.write("|")
        MyFile.write("\t")
        for a in B1[n]:
            MyFile.write(a)
            MyFile.write("\t")
        MyFile.write("\n")

    MyFile.write("\n\n")        
    MyFile.write("Counts of CLMotifs")        
    MyFile.write("\n")
    for n in B2:
        MyFile.write(n)
        MyFile.write("\t\t")
        MyFile.write("|")
        MyFile.write("\t")
        MyFile.write(str(B2[n]))
        MyFile.write("\n")
        
    MyFile.write("\n\n")         
    MyFile.write("hEXOmotifs")
    MyFile.write("\n")
    for n in C1:
        MyFile.write(n)
        MyFile.write("\t\t")
        MyFile.write("|")
        MyFile.write("\t")
        for a in C1[n]:
            MyFile.write(a)
            MyFile.write("\t")
        MyFile.write("\n")

    MyFile.write("\n\n")        
    MyFile.write("Counts of hEXOmotifs")        
    MyFile.write("\n")
    for n in C2:
        MyFile.write(n)
        MyFile.write("\t\t")
        MyFile.write("|")
        MyFile.write("\t")
        MyFile.write(str(C2[n]))
        MyFile.write("\n")
    MyFile.write("\n\n\n\n")         
    MyFile.close()

#Function: Main, compute all
def main():
    information_dict = opensequence(sys.argv[1])
    Exos=ExoMotif(information_dict)
    CLs=CLMotif(information_dict)
    hEXOs=hEXOmotifs(information_dict)
    savinginfo(Exos[0],Exos[1],CLs[0],CLs[1], hEXOs[0], hEXOs[1])                                                                                     
    

if __name__=="__main__":
    import sys#So as to use input arguments
    import re#So as to use regular expressions
    main()


