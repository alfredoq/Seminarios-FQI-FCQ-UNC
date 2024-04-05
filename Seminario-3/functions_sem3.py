import numpy as np
import pandas as pd
import sklearn
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem import rdDepictor
from IPython.display import SVG
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import IPythonConsole
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn import datasets
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
import seaborn as sns
import statsmodels.api as sm
import warnings
import py3Dmol
from itertools import combinations

def optimize_conf(mol):
    if mol is not None:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        return mol
    else:
        return None


def MolTo3DView(mol, size=(300, 300), style="stick", surface=False, opacity=0.5):
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
    mol_opt = optimize_conf(mol)
    mblock = Chem.MolToMolBlock(mol_opt)
    viewer = py3Dmol.view(width=size[0], height=size[1])
    viewer.addModel(mblock, 'mol')
    viewer.setStyle({style:{}})
    if surface:
        viewer.addSurface(py3Dmol.SAS, {'opacity': opacity})
    viewer.zoomTo()
    return viewer

    
def load_Smiles_supplier(smiles):
    suppl = Chem.SmilesMolSupplier(smiles, titleLine=False)
    return suppl


def molecules(template, suppl):
    template = Chem.MolFromSmiles(template)
    AllChem.Compute2DCoords(template)
    for mol in suppl:
        name = mol.GetProp("_Name")
        print(name)
        smile = Chem.MolToSmiles(mol)
        print(smile)
        AllChem.GenerateDepictionMatching2DStructure(mol, template)
        display(mol)
        viewer = MolTo3DView(mol, size=(300, 300), style="stick", surface=False, opacity=0.5)
        viewer.show()

def quiral_centers(suppl, template):
    template = Chem.MolFromSmiles(template)
    AllChem.Compute2DCoords(template)
    for mol in suppl:
        name = mol.GetProp("_Name")
        isomers = EnumerateStereoisomers(mol)
        LogP = Chem.Crippen.MolLogP(mol)
        for isomer in isomers:
            num_chiral_centers = Chem.FindMolChiralCenters(isomer, force=True)
            smile = Chem.MolToSmiles(isomer)
            AllChem.GenerateDepictionMatching2DStructure(isomer, template)
            display(isomer)
    
            print(smile)
            print(num_chiral_centers)
            print('LogP=',LogP)
            print(name)
    
def logP(mol):
    LogP = Chem.Crippen.MolLogP(mol)
    print('LogP=',LogP)

def rotable_bonds(mol):
    rotable_bonds = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)
    print('Enlaces_rotables', rotable_bonds)

def show_acid_base_groups(mol, size=(400, 200), kekulize=True):
    #mol = Chem.MolFromSmiles(smiles)
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

    #print(f"Basicos: {all_matches_basic}" )
    #print(f"Acidos: {all_matches_acid}" )
    all_matches = all_matches_basic + all_matches_acid

    
    drawer.DrawMolecule(mol, highlightAtoms=all_matches)
    
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return SVG(svg.replace('svg:',''))

def molecules_2(template, suppl):
    template = Chem.MolFromSmiles(template)
    AllChem.Compute2DCoords(template)
    for mol in suppl:
        name = mol.GetProp("_Name")
        print(name)
        smile = Chem.MolToSmiles(mol)
        print(smile)
        logP(mol)
        rotable_bonds(mol)
        show_acid_base_groups(mol)
        AllChem.GenerateDepictionMatching2DStructure(mol, template)
        display(mol)
        #viewer = MolTo3DView(mol, size=(300, 300), style="stick", surface=False, opacity=0.5)
        #viewer.show()