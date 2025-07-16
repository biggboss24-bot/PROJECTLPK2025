import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd

st.title('🪄 Reaction Pathway Predictor')
st.subheader('Prediksi Jalur Reaksi Kimia Berdasarkan Parameter Analisis')

# Sidebar untuk parameter input
with st.sidebar:
    st.header('Parameter Reaksi')
    temperature = st.slider('Suhu (°C)', 0, 200, 25)
    ph = st.slider('pH', 0, 14, 7)
    concentration = st.slider('Konsentrasi (mol/L)', 0.1, 5.0, 1.0)
    solvent = st.selectbox('Pelarut', ['Air', 'Etanol', 'Aseton', 'Dietil Eter', 'DMSO'])
    catalyst = st.checkbox('Ada Katalis?')
    
    if catalyst:
        catalyst_type = st.selectbox('Jenis Katalis', ['Asam', 'Basa', 'Logam Transisi', 'Enzim'])

# Input senyawa awal
st.header('Masukkan Reaktan')
col1, col2 = st.columns(2)
with col1:
    reactant1 = st.text_input('Reaktan 1 (SMILES)', 'C=O')
with col2:
    reactant2 = st.text_input('Reaktan 2 (SMILES)', 'O')

# Fungsi untuk validasi SMILES
def validate_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except:
        return False

# Fungsi untuk prediksi jalur reaksi (simulasi)
def predict_reaction_pathway(r1, r2, temp, ph, conc, solvent, catalyst=None):
    pathways = []
    if r1 == 'C=O' and r2 == 'O':
        pathways.append({
            'name': 'Hidrasi Aldehid',
            'steps': [
                ('Pembentukan intermediate tetrahedral', 50),
                ('Proton transfer', 30),
                ('Pembentukan produk hidrat', 20)
            ],
            'products': ['OC(O)C'],
            'activation_energy': 75 - (ph * 2) - (temp * 0.1),
            'thermodynamics': 'Eksotermik'
        })
    elif r1 == 'CC(=O)O' and r2 == 'CO':
        pathways.append({
            'name': 'Esterifikasi',
            'steps': [
                ('Protonasi karbonil', 60),
                ('Serangan nukleofilik', 40),
                ('Dehidrasi', 30)
            ],
            'products': ['CC(=O)OC', 'O'],
            'activation_energy': 90 - (ph * 3) - (temp * 0.2),
            'thermodynamics': 'Endotermik'
        })
    else:
        pathways.append({
            'name': 'Reaksi Adisi',
            'steps': [
                ('Inisiasi', 40),
                ('Propagasi', 30),
                ('Terminasi', 20)
            ],
            'products': [f'{r1}{r2}'],  # gabungan SMILES bisa jadi tidak valid!
            'activation_energy': 85 - (ph * 1.5) - (temp * 0.15),
            'thermodynamics': 'Eksotermik'
        })
    return pathways

# Tombol prediksi
if st.button('Prediksi Jalur Reaksi'):
    if not validate_smiles(reactant1) or not validate_smiles(reactant2):
        st.error('❌ Format SMILES tidak valid untuk satu atau lebih reaktan!')
    else:
        with st.spinner('🔍 Menganalisis kemungkinan jalur reaksi...'):
            try:
                pathways = predict_reaction_pathway(
                    reactant1, reactant2, temperature, ph, 
                    concentration, solvent, catalyst if catalyst else None
                )
                st.success('✅ Prediksi berhasil!')
                for pathway in pathways:
                    st.subheader(f"🔬 Jalur: {pathway['name']}")
                    
                    # Diagram energi
                    st.markdown("### ⚡ Diagram Energi Potensial")
                    fig, ax = plt.subplots()
                    steps = [step[0] for step in pathway['steps']]
                    energies = [step[1] for step in pathway['steps']]
                    steps = ['Reaktan'] + steps + ['Produk']
                    energies = [0] + energies + [energies[-1] - 20]
                    ax.plot(steps, energies, marker='o')
                    ax.set_ylabel('Energi (kJ/mol)')
                    ax.set_xlabel('Tahapan Reaksi')
                    ax.set_title('Profil Energi Reaksi')
                    plt.xticks(rotation=45)
                    plt.tight_layout()
                    st.pyplot(fig)

                    # Parameter kunci
                    col1, col2 = st.columns(2)
                    with col1:
                        st.metric("Energi Aktivasi", f"{pathway['activation_energy']:.1f} kJ/mol")
                    with col2:
                        st.metric("Termodinamika", pathway['thermodynamics'])

                    # Tahapan
                    st.markdown("### 🧪 Tahapan Reaksi")
                    for i, step in enumerate(pathway['steps'], 1):
                        st.markdown(f"{i}. *{step[0]}* - Energi: {step[1]} kJ/mol")

                    # Produk
                    st.markdown("### 🧫 Produk Prediksi")
                    mols = []
                    legends = []
                    for p in pathway['products']:
                        mol = Chem.MolFromSmiles(p)
                        if mol:
                            mols.append(mol)
                            legends.append(p)
                        else:
                            st.warning(f"⚠️ Produk tidak valid: {p}")

                    if mols:
                        img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(300, 300), legends=legends)
                        st.image(img)
                    else:
                        st.error("Tidak ada produk valid untuk divisualisasikan.")

                    st.divider()

            except Exception as e:
                st.error(f"❌ Terjadi kesalahan saat memproses: {e}")

# Analisis lanjutan
st.header('📊 Analisis Lanjutan')
with st.expander('🧠 Optimasi Kondisi Reaksi'):
    st.write("Rekomendasi berdasarkan parameter input:")
    data = {
        'Parameter': ['Suhu', 'pH', 'Konsentrasi', 'Pelarut'],
        'Nilai Saat Ini': [f"{temperature}°C", ph, f"{concentration} M", solvent],
        'Rekomendasi': [
            f"{temperature + 10}°C" if temperature < 100 else "Cukup",
            "Asam (pH 3-6)" if ph > 6 else "Basa (pH 8-11)" if ph < 8 else "Netral (baik)",
            f"{concentration * 1.2:.1f} M" if concentration < 3 else "Cukup",
            "Pertahankan" if solvent in ['Air', 'Etanol'] else "Coba ganti ke Air"
        ]
    }
    st.table(pd.DataFrame(data))

with st.expander('📄 Simpan Hasil Analisis'):
    st.download_button(
        label="📥 Unduh Laporan PDF",
        data="Simulasi laporan PDF. Akan berisi data aktual dalam versi produksi.",
        file_name="reaction_analysis.pdf",
        mime="application/pdf"
    )

# Footer
st.caption("""
🧪 Aplikasi ini menggunakan model prediksi sederhana untuk keperluan demonstrasi. 
Hasil aktual bisa berbeda bergantung pada kondisi eksperimen nyata.
""")
