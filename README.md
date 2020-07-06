# generate_initial
Created by Tim Ling @ YSL lab

## Dependencies

### Python
* numpy
* scipy

### MDAnalysis
* MDAnalysis

### Chimera
* Chimera 1.14 is used when testing this program.

## Protocol Description
	
	This program generates initial structures for MD simulations.

## Quickstart

  To run the program on the cluster, ensure that MDAnalysis is installed and Chimera is loaded.
  You will need to change 
  
  	"/Applications/Chimera.app/Contents/MacOS/chimera" 
	
  in chimScriptMaker.py to the binary of your installation of simply chimera if the module is loaded properly. 
	
  To run the program, simply do:

				python main.py --seq SEQUENCE 

	The program can also take the following flags:
		
    --cutoffRMSD: (float) the cutoff for RMSD between all structures. Default is 1.0.
    
    --cutoffOmega: (float) the cutoff for Omega angle. Default is 150, which means a peptide bond is considered 
                   as trans if the omega angle has an absolute value greaer than 150 and cis if otherwise.
    
    --cutoffBond: (float) the cutoff for head-to-tail distance. Default is 1.4 Angstrom.
    
    --n: (int) number of initial structures to be generated. Default is 2
    
    --cyclic: (bool) if the structure is cyclic. 
    
## Important Note:

  The program is currently still under testing phase and the omega angle, chirality, and head-to-tail distance should be checked manually.
  
## Examples:
### SESEaaDG-s1 with Omega angles labelled
![](https://github.com/yling01/generate_initial/blob/master/Picture/SESEaaDG.png)


## Acknowledgments
This code is based on the first version written by Aidan Fike and was later updated by Kevin Schult.




