from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDepictor
from IPython.display import SVG
import py3Dmol
import numpy as np
import matplotlib.pyplot as plt


def construct_molecule_2d(smiles):
    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol)
    display(mol)

    return mol

def optimize_conf(mol):
    if mol is not None:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        return mol
    else:
        return None

def MolTo3DView(smiles, size=(300, 300), style="stick", surface=False, opacity=0.5):
    """Draw molecule in 3D
    
    Args:
    ----
        mol: rdMol, molecule to show
        size: tuple(int, int), canvas size
        style: str, type of drawing molecule
               style can be 'line', 'stick', 'sphere', 'carton'
        surface, bool, display SAS
        opacity, float, opacity of surface, range 0.0-1.0
    Return:
    ----
        viewer: py3Dmol.view, a class for constructing embedded 3Dmol.js views in ipython notebooks.
    """
    assert style in ('line', 'stick', 'sphere', 'carton')
    mol = Chem.MolFromSmiles(smiles)
    mol_opt = optimize_conf(mol)
    mblock = Chem.MolToMolBlock(mol_opt)
    viewer = py3Dmol.view(width=size[0], height=size[1])
    viewer.addModel(mblock, 'mol')
    viewer.setStyle({style:{}})
    if surface:
        viewer.addSurface(py3Dmol.SAS, {'opacity': opacity})
    viewer.zoomTo()
    return viewer

def show_acid_base_groups(smiles, size=(400, 200), kekulize=True):
    mol = Chem.MolFromSmiles(smiles)
    rdDepictor.Compute2DCoords(mol)
    kekulize = True
    if kekulize:
        Chem.Kekulize(mol) # Localize the benzene ring bonds
            
    drawer = rdMolDraw2D.MolDraw2DSVG(400, 200)

    basic_groups = {'primary_amine' : '[NX3;H2;!$(NC=O)]', # This excludes amides
                'secondary_amine' :'[NX3;H1;!$(NC=O)]', # This excludes amides
                'tertiary_amine' : '[NX3;H0;!$(NC=O)]'} # This excludes amides 
      
    acid_groups = {'sulphonamide' : '[SX4]([NX3])(=[OX1])(=[OX1])',
               'carboxylic_acid' : '[CX3](=O)[OX1H0-,OX2H1]',
               'phenolic_hydroxyl' : '[$([OX2;H1][cX3])]',
               'barbituric_acid' : '[CX3](=[O])[NX3;H1][CX3](=[O])',}
    
    all_matches_basic = ()
    all_matches_acid = ()

    for key in basic_groups:
        group = basic_groups[key]
        substructure = Chem.MolFromSmarts(group)
        # highlightAtoms expects only one tuple, not tuple of tuples. So it needs to be merged into a single tuple
        matches = sum(mol.GetSubstructMatches(substructure), ())
        #print(key, matches)
        all_matches_basic += matches

    for key in acid_groups:
        group = acid_groups[key]
        substructure = Chem.MolFromSmarts(group)
        # highlightAtoms expects only one tuple, not tuple of tuples. So it needs to be merged into a single tuple
        matches = sum(mol.GetSubstructMatches(substructure), ())
        #print(key, matches)
        all_matches_acid += matches

    print(f"Basicos: {all_matches_basic}" )
    print(f"Acidos: {all_matches_acid}" )
    all_matches = all_matches_basic + all_matches_acid

    
    drawer.DrawMolecule(mol, highlightAtoms=all_matches)
    
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return SVG(svg.replace('svg:',''))


def calculate_ionized_fraction(pka,specie='acid',ph1=1, ph2=14,step=0.5,return_data=0):
    ph_values_list = []
    ionized_fraction_list = []
    neutral_fraction_list = []
    
    for ph in np.arange(ph1, ph2+step, step):
        if specie == 'acid':
            ionized_fraction = round((np.exp(ph - pka) * 100) / (np.exp(ph - pka) + 1),2)
            neutral_fraction = round(100-ionized_fraction,2)
            if return_data == 0:
                print(f"pH: {ph}")
                print(f"Especie: {specie} - porcentaje de base conjugada ionizada [A-]: {ionized_fraction} %")
                print(f"Especie: {specie} - porcentaje de especie neutra [AH-]: {neutral_fraction} % \n")
                print(f'-------------------------------------------------------------------------')  
        elif specie == 'base':
            ionized_fraction = round(100 - ((np.exp(ph - pka) * 100) / (np.exp(ph - pka) + 1)),2)
            neutral_fraction = round(100-ionized_fraction,2)
            if return_data == 0:
                print(f"pH: {ph}")
                print(f"Especie: {specie} - porcentaje de ácido conjugado ionizada [BH+]: {ionized_fraction} %")
                print(f"Especie: {specie} - porcentaje de especie neutra [B]: {neutral_fraction} % \n")
                print(f'-------------------------------------------------------------------------')
        else:
            print("Indicar si es 'acid' o 'base'")
            return

        ph_values_list.append(ph)
        ionized_fraction_list.append(ionized_fraction)
        neutral_fraction_list.append(neutral_fraction)

    if return_data == 1:
        return ph_values_list, ionized_fraction_list, neutral_fraction_list


def plot_ionized_fraction(pka,specie='acid',ph1=1, ph2=14,step=0.5):
    if specie == 'acid':
        ph_values_list, ionized_fraction_list, neutral_fraction_list = calculate_ionized_fraction(pka,specie,ph1, ph2,step,return_data=1)
        plt.scatter(ph_values_list,ionized_fraction_list)
        plt.plot(ph_values_list,ionized_fraction_list,c='r')
        plt.xticks(np.arange(ph1, ph2+step, step),rotation='vertical')
        plt.yticks(np.arange(0, 100, 10))
        plt.title("Ionized fraction of [A-] for specie ACID")
        plt.xlabel("pH")
        plt.ylabel("[A-] %")

    elif specie == 'base':
        ph_values_list, ionized_fraction_list, neutral_fraction_list = calculate_ionized_fraction(pka,specie,ph1, ph2,step,return_data=1)
        
        # Plot the ionized fraction
        plt.scatter(ph_values_list,ionized_fraction_list,label='Ionizada',c='r')
        plt.plot(ph_values_list,ionized_fraction_list,c='r')
        
        # Plot the neutral fraction
        plt.scatter(ph_values_list,neutral_fraction_list,label='Neutra',c='g',alpha=0.5)
        plt.plot(ph_values_list,neutral_fraction_list,c='g',alpha=0.5)
                
        plt.xticks(np.arange(ph1, ph2+step, step),rotation='vertical')
        plt.yticks(np.arange(0, 100, 10))
        plt.title("Ionized fraction of [BH+] for specie BASE")
        plt.xlabel("pH")
        plt.ylabel("[BH+] %")
        plt.legend()
        
    plt.grid()
    plt.show()
    


    
    