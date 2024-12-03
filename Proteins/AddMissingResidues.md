#Adding missing residues in files downloaded from PDB databank (simulation preparation)


When you download a PDB, READ the HEADER. You find a number of useful information including (i) the number of identical chains (we typically do simulations on one chain) (ii) missing residues. For example, if you download the PDB file ID:6W63, you will find the following remark:

```
REMARK 465 MISSING RESIDUES                                                     
REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE                       
REMARK 465 EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN               
REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)                
REMARK 465                                                                      
REMARK 465   M RES C SSSEQI                                                     
REMARK 465     GLN A   306                                                      
REMARK 500                                                                      
```
meaning that residue #306 is missing. You can double check this information by looking at the sequence of the protein (also in the header):

```
SEQRES   1 A  306  SER GLY PHE ARG LYS MET ALA PHE PRO SER GLY LYS VAL          
SEQRES   2 A  306  GLU GLY CYS MET VAL GLN VAL THR CYS GLY THR THR THR          
SEQRES   3 A  306  LEU ASN GLY LEU TRP LEU ASP ASP VAL VAL TYR CYS PRO          
SEQRES   4 A  306  ARG HIS VAL ILE CYS THR SER GLU ASP MET LEU ASN PRO          
SEQRES   5 A  306  ASN TYR GLU ASP LEU LEU ILE ARG LYS SER ASN HIS ASN          
SEQRES   6 A  306  PHE LEU VAL GLN ALA GLY ASN VAL GLN LEU ARG VAL ILE          
SEQRES   7 A  306  GLY HIS SER MET GLN ASN CYS VAL LEU LYS LEU LYS VAL          
SEQRES   8 A  306  ASP THR ALA ASN PRO LYS THR PRO LYS TYR LYS PHE VAL          
SEQRES   9 A  306  ARG ILE GLN PRO GLY GLN THR PHE SER VAL LEU ALA CYS          
SEQRES  10 A  306  TYR ASN GLY SER PRO SER GLY VAL TYR GLN CYS ALA MET          
SEQRES  11 A  306  ARG PRO ASN PHE THR ILE LYS GLY SER PHE LEU ASN GLY          
SEQRES  12 A  306  SER CYS GLY SER VAL GLY PHE ASN ILE ASP TYR ASP CYS          
SEQRES  13 A  306  VAL SER PHE CYS TYR MET HIS HIS MET GLU LEU PRO THR          
SEQRES  14 A  306  GLY VAL HIS ALA GLY THR ASP LEU GLU GLY ASN PHE TYR          
SEQRES  15 A  306  GLY PRO PHE VAL ASP ARG GLN THR ALA GLN ALA ALA GLY          
SEQRES  16 A  306  THR ASP THR THR ILE THR VAL ASN VAL LEU ALA TRP LEU          
SEQRES  17 A  306  TYR ALA ALA VAL ILE ASN GLY ASP ARG TRP PHE LEU ASN          
SEQRES  18 A  306  ARG PHE THR THR THR LEU ASN ASP PHE ASN LEU VAL ALA          
SEQRES  19 A  306  MET LYS TYR ASN TYR GLU PRO LEU THR GLN ASP HIS VAL          
SEQRES  20 A  306  ASP ILE LEU GLY PRO LEU SER ALA GLN THR GLY ILE ALA          
SEQRES  21 A  306  VAL LEU ASP MET CYS ALA SER LEU LYS GLU LEU LEU GLN          
SEQRES  22 A  306  ASN GLY MET ASN GLY ARG THR ILE LEU GLY SER ALA LEU          
SEQRES  23 A  306  LEU GLU ASP GLU PHE THR PRO PHE ASP VAL VAL ARG GLN          
SEQRES  24 A  306  CYS SER GLY VAL THR PHE GLN 

```
where residue #305 is PHE and residue #306 is GLN. However, in the atom lines, we have:

```sh
ATOM   4637  CG2 THR A 304      13.740  31.548 -35.016  1.00 24.87           C  
ANISOU 4637  CG2 THR A 304     2941   3640   2869   -394   -143    295       C  
ATOM   4638  N   PHE A 305      11.263  28.467 -36.891  1.00 50.82           N  
ANISOU 4638  N   PHE A 305     6225   6934   6151   -319    -42    171       N  
ATOM   4639  CA  PHE A 305      10.653  27.515 -37.838  1.00 59.54           C  
ANISOU 4639  CA  PHE A 305     7321   8061   7242   -306    -16    136       C  
ATOM   4640  C   PHE A 305      10.370  28.183 -39.185  1.00 60.66           C  
ANISOU 4640  C   PHE A 305     7447   8247   7353   -357    -38    163       C  
ATOM   4641  O   PHE A 305      10.992  27.851 -40.203  1.00 55.63           O  
ANISOU 4641  O   PHE A 305     6774   7690   6672   -373    -25    159       O  
ATOM   4642  CB  PHE A 305       9.361  26.916 -37.247  1.00 55.24           C  
ANISOU 4642  CB  PHE A 305     6810   7439   6740   -268     -5    105       C  
ATOM   4643  CG  PHE A 305       9.492  26.545 -35.797  1.00 56.82           C  
ANISOU 4643  CG  PHE A 305     7029   7588   6972   -229      6     91       C  
ATOM   4644  CD1 PHE A 305       9.953  25.294 -35.427  1.00 54.28           C  
ANISOU 4644  CD1 PHE A 305     6700   7273   6652   -188     39     54       C  
ATOM   4645  CD2 PHE A 305       9.205  27.474 -34.800  1.00 54.86           C  
ANISOU 4645  CD2 PHE A 305     6804   7287   6752   -232    -20    114       C  
ATOM   4646  CE1 PHE A 305      10.105  24.961 -34.092  1.00 54.15           C  
ANISOU 4646  CE1 PHE A 305     6700   7211   6663   -155     47     46       C  
ATOM   4647  CE2 PHE A 305       9.361  27.153 -33.462  1.00 49.82           C  
ANISOU 4647  CE2 PHE A 305     6181   6608   6138   -198     -9    102       C  
ATOM   4648  CZ  PHE A 305       9.819  25.899 -33.107  1.00 53.30           C  
ANISOU 4648  CZ  PHE A 305     6614   7057   6579   -162     24     71       C  
TER    4649      PHE A 305                                                      
HETATM 4650  C02 X77 A 401     -19.908  21.124 -29.295  1.00 31.51           C 
``` 

where the last residue is PHE 305, followed by the complexed inhibitor (X77). Residue #306 is indeed missing. 

If there are multiple chains, choose the one with the least missing residues or the best ligand geometry (do check all chains in VMD or similar to see if there are obvious differences between the ligands/inhibitor).

Here, we will continue with chain A of 6W63. Since we know the sequence, we can add the missing residue to the file. We do that with MODELLER. 


###Install Modeller 

Install MODELLER on your laptop following the instructions [here](https://salilab.org/modeller/download_installation.html). It is free to academics but you will have to register for a license. 

### Write the sequence to an alignment file

Create a directory on your local machine where you will do this work and copy the pdb you have downloaded from the databank. 

```sh
cd wherever_you_want_on_your_laptop
mkdir my_new_directory
cd my_new_directory
cp ~/Downloads/6w63.pdb ./
```

Add the Python script `GetSequence.py` to `my_new_directory`:

```sh
from modeller import *
code = '6w63'

e = environ()
m = model(e, file=code)
aln = alignment(e)
aln.append_model(m, align_codes=code)
aln.write(file=code+'.seq')

```
Execute with

```sh
python GetSequence.py
```

This will write a sequence file, `6w63.seq` that contains the following:

```sh

>P1;6w63
structureX:6w63:   1 :A:+305 :A:MOL_ID  1; MOLECULE  MAIN PROTEASE; CHAIN  A; ENGINEERED  YES:MOL_ID  1; ORGANISM_SCIENTIFIC  SARS-COV-2; ORGANISM_TAXID  2697049; STRAIN  COVID-19; EXPRESSION_SYSTEM  ESCHERICHIA COLI BL21(DE3); EXPRESSION_SYSTEM_TAXID  469008; EXPRESSION_SYSTEM_STRAIN  BL21(DE3): 2.10: 0.16
SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDVVYCPRHVICTSEDMLNPNYEDLLIRKSNHNFLVQAGNVQL
RVIGHSMQNCVLKLKVDTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIKGSFLNGSCGSVGF
NIDYDCVSFCYMHHMELPTGVHAGTDLEGNFYGPFVDRQTAQAAGTDTTITVNVLAWLYAAVINGDRWFLNRFTT
TLNDFNLVAMKYNYEPLTQDHVDILGPLSAQTGIAVLDMCASLKELLQNGMNGRTILGSALLEDEFTPFDVVRQC
SGVTF*
```

From this file, create the file `alignment.ali` that contains a copy of the above sequence with `-` in place of missing residues and a copy of the full sequence. In our example, only the last residue is missing so we get: (notice the change of name in the header above the second sequence). 

```sh
>P1;6w63
structureX:6w63:   1 :A:+305 :A:MOL_ID  1; MOLECULE  MAIN PROTEASE; CHAIN  A; ENGINEERED  YES:MOL_ID  1; ORGANISM_SCIENTIFIC  SARS-COV-2; ORGANISM_TAXID  2697049; STRAIN  COVID-19; EXPRESSION_SYSTEM  ESCHERICHIA COLI BL21(DE3); EXPRESSION_SYSTEM_TAXID  469008; EXPRESSION_SYSTEM_STRAIN  BL21(DE3): 2.10: 0.16
SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDVVYCPRHVICTSEDMLNPNYEDLLIRKSNHNFLVQAGNVQL
RVIGHSMQNCVLKLKVDTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIKGSFLNGSCGSVGF
NIDYDCVSFCYMHHMELPTGVHAGTDLEGNFYGPFVDRQTAQAAGTDTTITVNVLAWLYAAVINGDRWFLNRFTT
TLNDFNLVAMKYNYEPLTQDHVDILGPLSAQTGIAVLDMCASLKELLQNGMNGRTILGSALLEDEFTPFDVVRQC
SGVTF-*
>P1;6w63_fill
sequence:::::::::
SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDVVYCPRHVICTSEDMLNPNYEDLLIRKSNHNFLVQAGNVQL
RVIGHSMQNCVLKLKVDTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIKGSFLNGSCGSVGF
NIDYDCVSFCYMHHMELPTGVHAGTDLEGNFYGPFVDRQTAQAAGTDTTITVNVLAWLYAAVINGDRWFLNRFTT
TLNDFNLVAMKYNYEPLTQDHVDILGPLSAQTGIAVLDMCASLKELLQNGMNGRTILGSALLEDEFTPFDVVRQC
SGVTFQ*

```


###Add the missing residue
We use the automodel routine of MODELLER which, by default re-optimize all the residues after adding the missing one(s). Since we have a complexed inhibitor here (not included in MODELLER algorithms), we will keep all residues fixed and simply add one at the end to keep the docking geometry from the orginal PDB file. 

`AddMissing.py`:

```sh
from modeller import *
from modeller.automodel import *    # Load the automodel class
   
log.verbose()
env = environ()
    
# directories for input atom files
env.io.atom_files_directory = ['.']
    
class MyModel(automodel):
    def select_atoms(self):
          return selection(self.residue_range('306', '306'))
   
a = MyModel(env, alnfile = 'alignment.ali',
               knowns = '6w63', sequence = '6w63_fill')
a.starting_model= 1
a.ending_model  = 1
   
a.make()
```


Execute with

```sh
python AddMissing.py
```

The resulting model structure is `6w63_fill.B99990001.pdb`. It is worth checking in VMD or similar visualization software to make sure you are happy with the changes. 

If you had a ligand (like here inhibitor X77), grep the `X77` lines and cat them the end of 6w63_fill.B99990001.pdb to get a complete PDB.

Note that MODELLER removed structural waters and protons. Structural waters shouldn't be added again but we will add protons.