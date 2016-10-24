package main

import (
	"fmt"
	"github.com/rmera/gochem"
	"os"
	"strings"
)


/*This program takes a PDB protonated with the Reduce program. It fixes duplicated atoms by leaving only
the copy labeled "A" in the crystallographic structure (assumed to be the one with the highest occupancy)
and records the protonation state of each ASP, GLU, LYS, ARG and HIS residue. These protonation states are
then printed as a string of "\n"-separated numbers, so taht string can be fed into the gromacs program 
gmx pdb2gmx ran with the options -his -arg -lys -asp -glu to ensure that the resulting file has the same protonation
states as those obtained with Reduce. The program will also output a PDB without the hydrogens, which can be
used as an input for gmx pdb2gmx*/


//If there are four atoms with names containing "HH", the ARG is considered protonated
func Arg(mol *chem.Molecule, resi, i int) (newi, protonated int) {
	protonated = 0
	HH := 0 //if it gets to four, the ARG is protonated
	for newi = i; newi < mol.Len(); newi++ {
		at := mol.Atom(newi)
		if at.MolID != resi {
			break
		}
		if strings.HasPrefix(at.Name, "HH") {
			HH++
		}
	}
	if HH == 4 {
		protonated = 1
	}
	newi--
	return

}

//if there are 3 atoms with names containing "HZ", the ARG is considered protonated
func Lys(mol *chem.Molecule, resi, i int) (newi, protonated int) {
	protonated = 0
	HZ := 0 //if it gets to four, the ARG is protonated
	for newi = i; newi < mol.Len(); newi++ {
		at := mol.Atom(newi)
		if at.MolID != resi {
			break
		}
		if strings.HasPrefix(at.Name, "HZ") {
			HZ++
		}
	}
	if HZ == 3 {
		protonated = 1
	}
	newi--
	return

}

//1 if protonated, 0 if not
func GluAsp(mol *chem.Molecule, resi, i int) (newi, protonated int) {
	protonated = 0
	comp := "HD" //The string to compare to the atoms name, "HD" if ASP, HE if GLU. I am not 100% sure that these names are correct.
	if mol.Atom(i).Molname == "GLU" {
		comp = "HE"
	}
	for newi = i; newi < mol.Len(); newi++ {
		at := mol.Atom(newi)
		if at.MolID != resi {
			newi--
			return
		}
		if strings.HasPrefix(at.Name, comp) {
			fmt.Println("Name in", at, at.Name, "Coords", mol.Coords[0].VecView(newi))
			protonated = 1
			newi--
			return
		}
	}
	panic("unreachable")
}

func His(mol *chem.Molecule, resi, i int) (newi, protonated int) {
	protonated = 0
	HD1 := false //if it gets to four, the ARG is protonated
	HE2 := false
	for newi = i; newi < mol.Len(); newi++ {
		at := mol.Atom(newi)
		if at.MolID != resi {
			break
		}
		if strings.HasPrefix(at.Name, "HD1") {
			HD1 = true
		}
		if strings.HasPrefix(at.Name, "HE2") {
			HE2 = true
		}
	}
	if HD1 && HE2 {
		protonated = 2
	} else if HE2 {
		protonated = 1
	} else if HD1 {
		protonated = 0
	} else {
		panic("Abnormal histidine! perhaps Heme-coupled or as Imidazolate (SOD1)?")
	}
	newi--
	return

}


func main() {
	mol, err := chem.PDBFileRead(os.Args[3], false)
	if err != nil {
		panic("Error reading PDB file: " + err.Error())
	}
	//we'll get the chains also
	Chains := make([]string, 0, 4)
	//We'll start by cleaning up duplicated atoms, we simply stick with the first altenative
	for i := 0; i < mol.Len(); i++ {
		at := mol.Atom(i)
		if i == 0 {
			Chains = append(Chains, at.Chain)
		} else if at.Chain != Chains[len(Chains)-1] {
			Chains = append(Chains, at.Chain)
		}
		//	fmt.Println(string(at.Char16))
		//We keep the first version of the atom (A) and delete the others. I assumes there are no more than 4 versions! the "C" could interfere if used for a terminal residue.
		//This will also fail if non-unicode names are used for atoms.
		if at.Char16 == 'B' || at.Char16 == 'C' || at.Char16 == 'D' {
			//		fmt.Println(at, at.Char16)
			mol.Del(i)
			i--
		}
	}
	asp := make([][2]int, 0, 10)
	glu := make([][2]int, 0, 10)
	his := make([][2]int, 0, 10)
	lys := make([][2]int, 0, 10)
	arg := make([][2]int, 0, 10)
	var ResidueList *[][2]int
	//Now we create the input for pdb2gmx -his -asp -glu -arg -lys
	prevres := -1
	gmx := fmt.Sprintf("%s\n%s\n", os.Args[1], os.Args[2]) //first is the number for the wanted force-field, second is the number for the wanted water model.
	var Protonator func(*chem.Molecule, int, int) (int, int)
	for i := 0; i < mol.Len(); i++ {
		at := mol.Atom(i)
		if at.MolID == prevres {
			continue
		}
		curres := at.MolID
		prevres = at.MolID
		prot := 0
		switch at.Molname {
		case "HIS":
			Protonator = His
			ResidueList = &his
		case "GLU":
			Protonator = GluAsp
			ResidueList = &glu
		case "ASP":
			Protonator = GluAsp
			ResidueList = &asp
		case "ARG":
			Protonator = Arg
			ResidueList = &arg
		case "LYS":
			Protonator = Lys
			ResidueList = &lys
		default:
			Protonator = nil
			ResidueList = nil

		}
		if Protonator == nil {
			continue
		}
		i, prot = Protonator(mol, curres, i)
		fmt.Println([2]int{i, prot})
		*ResidueList = append(*ResidueList, [2]int{i, prot}) 	
	}
	for _, chain := range Chains {
		for _, v := range [][][2]int{lys, arg, asp, glu, his} {  //This residue order is the one Gromacs uses
			for _, v2 := range v {
				at := mol.Atom(v2[0])
				if at.Chain != chain {
					continue
				}
				fmt.Fprintf(os.Stderr, "Residue: %s %d, of Chain %s, assigned gmx protonation code %d\n", at.Molname, at.MolID, at.Chain, v2[1])
				gmx = fmt.Sprintf("%s%d\n", gmx, v2[1])
			}
		}
	}

	chem.PDBFileWrite(strings.Replace(os.Args[3], ".pdb", "_clean.pdb", 1), mol.Coords[0], mol, nil)

	for i := 0; i < mol.Len(); i++ {
		at := mol.Atom(i)
		if at.Symbol == "H" {
			mol.Del(i)
			i--
		}
	}

	chem.PDBFileWrite(strings.Replace(os.Args[3], ".pdb", "_cleanNoH.pdb", 1), mol.Coords[0], mol, nil)

	fmt.Printf("%s", gmx)
	fmt.Println(lys, arg, glu, his, asp)
}
