import numpy as np
import pandas as pd
import streamlit as st
import pickle
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Descriptors

st.write("# Molecular Solubility Predictor")
st.write("")

image = Image.open('solubility-logo.jpg')

st.image(image, use_column_width=True)

st.write("")
st.write("""

This app predicts the **Solubility (LogS)** values of molecules!

Data obtained from the John S. Delaney. [ESOL:  Estimating Aqueous Solubility Directly from Molecular Structure](https://pubs.acs.org/doi/10.1021/ci034243x). ***J. Chem. Inf. Comput. Sci.*** 2004, 44, 3, 1000-1005.
***
""")

## Calculate molecular descriptors
def AromaticProportion(m):
  aromatic_atoms = [m.GetAtomWithIdx(i).GetIsAromatic() for i in range(m.GetNumAtoms())]
  aa_count = []
  for i in aromatic_atoms:
    if i==True:
      aa_count.append(1)
  AromaticAtom = sum(aa_count)
  HeavyAtom = Descriptors.HeavyAtomCount(m)
  AR = AromaticAtom/HeavyAtom
  return AR

def generate(smiles, verbose=False):

    moldata= []
    for elem in smiles:
        mol=Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData= np.arange(1,1)
    i=0
    for mol in moldata:

        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_MolWt = Descriptors.MolWt(mol)
        desc_NumRotatableBonds = Descriptors.NumRotatableBonds(mol)
        desc_AromaticProportion = AromaticProportion(mol)

        row = np.array([desc_MolLogP,
                        desc_MolWt,
                        desc_NumRotatableBonds,
                        desc_AromaticProportion])

        if(i==0):
            baseData=row
        else:
            baseData=np.vstack([baseData, row])
        i=i+1

    columnNames=["MolLogP","MolWt","NumRotatableBonds","AromaticProportion"]
    descriptors = pd.DataFrame(data=baseData,columns=columnNames)

    return descriptors

st.sidebar.title('Input Canonical Smiles')
st.sidebar.write('Molcular Solubility Predictor, a comprehensive tool to predict solubility of discovered as well as undiscovered molecules.')
st.write("")

## Read SMILES input
sample_input = "NCCCC\nCCC\nCN"

SMILES = st.sidebar.text_area("User Input", sample_input)

def update(SMILES):
    SMILES = "C\n" + SMILES #Adds C as a dummy, first item
    SMILES = SMILES.split('\n')

    st.header('Input SMILES')
    SMILES[1:] # Skips the dummy first item

    ## Calculate molecular descriptors
    st.header('Computed molecular descriptors')
    X = generate(SMILES)
    X[1:] # Skips the dummy first item

    # Reads in saved model
    load_model = pickle.load(open('solubility_model.pkl', 'rb'))

    # Apply model to make predictions
    prediction = load_model.predict(X)

    st.header('Predicted LogS values')
    predictions_with_smiles = pd.DataFrame({'SMILES': SMILES[1:], 'Predicted LogS': prediction[1:],})
    st.write(predictions_with_smiles)

predict = st.sidebar.button('Predict', on_click=update(SMILES))