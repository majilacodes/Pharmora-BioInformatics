import pandas as pd

sol = pd.read_csv('delaney.csv')
sol

from rdkit import Chem

mol_list = [Chem.MolFromSmiles(element) for element in sol.SMILES]

len(mol_list)

mol_list[:5]

import numpy as np
from rdkit.Chem import Descriptors

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

X = generate(sol.SMILES)

X

Y = sol.iloc[:,1]
Y = Y.rename("logS")
Y

Y

Y.hist()

dataset = pd.concat([X,Y], axis=1)
dataset

dataset.to_csv('delaney_solubility_with_descriptors.csv', index=False)

X = dataset.drop(['logS'], axis=1)
X

Y = dataset.iloc[:,-1]
Y

"""#Linear Regression Model"""

from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score

model = linear_model.LinearRegression()
model.fit(X, Y)

Y_pred = model.predict(X)
Y_pred

# Commented out IPython magic to ensure Python compatibility.
a = mean_squared_error(Y, Y_pred)
b = r2_score(Y, Y_pred)

print(f'Coefficients: {model.coef_}')
print(f'Intercept: {model.intercept_}')

print(f'Mean squared error (MSE): {a:.2f}')

print(f'Coefficient of determination (R^2): {b:.2f}')

"""#Model Equation"""

print('LogS = %.2f %.2f LogP %.4f MW + %.4f RB %.2f AP' % (model.intercept_, model.coef_[0], model.coef_[1], model.coef_[2], model.coef_[3] ) )

import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=(5,5))
plt.scatter(x=Y, y=Y_pred, c="#7CAE00", alpha=0.3)

# Add trendline
# https://stackoverflow.com/questions/26447191/how-to-add-trendline-in-python-matplotlib-dot-scatter-graphs
z = np.polyfit(Y, Y_pred, 1)
p = np.poly1d(z)

plt.plot(Y,p(Y),"#F8766D")
plt.ylabel('Predicted LogS')
plt.xlabel('Experimental LogS')
plt.savefig('scatter_plot.pdf', format='pdf')

import pickle

pickle.dump(model, open('solubility_model.pkl', 'wb'))


