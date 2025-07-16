import streamlit as st
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd

st.title('🪄 Reaction Pathway Predictor')
st.subheader('Prediksi Jalur Reaksi Kimia Berdasarkan Parameter Analisis')

# Sidebar parameter
with st.sidebar:
    st.header('Parameter Reaksi')
    temperature = st.slider('Suhu (°C)', 0, 200, 25)
    ph = st.slider('pH', 0, 14, 7)
    concentration = st.slider('Konsentrasi (mol/L)', 0.1, 5.0, 1.0)
    solvent = st.selectbox('Pelarut', ['Air', 'Etanol', 'Aseton', 'Dietil Eter', 'DMSO'])
    catalyst = st.checkbox('Ada Katalis?')
    if catalyst:
        catalyst_type = st.selectbox('Jenis Katalis', ['Asam', 'Basa', 'Logam Transisi', 'Enzim'])

# Input reaktan
st.header('Masukkan Reaktan')
col1, col2 = st.columns(2)
with col1:
    reactant1 = st.text_input('Reaktan 1 (SMILES)', 'C=O')
with col2:
    reactant2 = st.text_input('Reaktan 2 (SMILES)', 'O')

# Validasi SMILES
def validate_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except:
        return False

# Prediksi jalur reaksi (simulasi)
def predict_reaction_pathway(r1, r2, temp, ph, conc, solvent, catalyst=None):
    pathways = []
    if r1 == 'C=O' and r2 == 'O':
        pathways.append({
            'name': 'Hidrasi Aldehid',
            'steps': [
                ('Intermediate tetrahedral', 50),
                ('Proton transfer', 30),
                ('Hidrat terbentuk'
