import streamlit as st
import pandas as pd

st.title('🪄 Reaction Pathway Predictor (Tanpa RDKit)')
st.subheader('Prediksi Jalur Reaksi Kimia Berdasarkan Parameter Analisis')

# Sidebar input parameter
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

# Prediksi jalur reaksi
def predict_reaction_pathway(r1, r2, temp, ph, conc, solvent, catalyst=None):
    pathways = []
    if r1 == 'C=O' and r2 == 'O':
        pathways.append({
            'name': 'Hidrasi Aldehid',
            'steps': [
                ('Intermediate tetrahedral terbentuk', 50),
                ('Proton transfer', 30),
                ('Produk hidrat terbentuk', 20)
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
            'products': ['CC(=O)OC + H2O'],
            'activation_energy': 90 - (ph * 3) - (temp * 0.2),
            'thermodynamics': 'Endotermik'
        })
    else:
        pathways.append({
            'name': 'Reaksi Adisi Umum',
            'steps': [
                ('Inisiasi', 40),
                ('Propagasi', 30),
                ('Terminasi', 20)
            ],
            'products': [f'{r1}.{r2}'],
            'activation_energy': 85 - (ph * 1.5) - (temp * 0.15),
            'thermodynamics': 'Eksotermik'
        })
    return pathways

# Tombol prediksi
if st.button('Prediksi Jalur Reaksi'):
    if not reactant1.strip() or not reactant2.strip():
        st.error('❌ Mohon masukkan dua reaktan dalam format SMILES.')
    else:
        with st.spinner('🔬 Menganalisis...'):
            results = predict_reaction_pathway(
                reactant1, reactant2, temperature, ph, concentration, solvent, catalyst
            )
            st.success('✅ Prediksi berhasil!')
            for path in results:
                st.subheader(f"🔹 Jalur Reaksi: {path['name']}")

                st.markdown("**Energi Aktivasi:** {:.1f} kJ/mol".format(path['activation_energy']))
                st.markdown("**Termodinamika:** {}".format(path['thermodynamics']))

                st.markdown("### 🧪 Tahapan Reaksi")
                for i, (step, energy) in enumerate(path['steps'], 1):
                    st.write(f"{i}. {step} — {energy} kJ/mol")

                st.markdown("### 🧫 Produk Prediksi")
                for prod in path['products']:
                    st.code(prod)

                st.divider()

# Optimasi kondisi
st.header("📊 Optimasi Reaksi")
with st.expander("Rekomendasi Kondisi"):
    data = {
        'Parameter': ['Suhu', 'pH', 'Konsentrasi', 'Pelarut'],
        'Nilai Saat Ini': [f"{temperature}°C", ph, f"{concentration} M", solvent],
        'Rekomendasi': [
            f"{temperature + 10}°C" if temperature < 100 else "Cukup",
            "Gunakan kondisi asam (pH 3-6)" if ph > 7 else "Coba pH basa (8-11)" if ph < 6 else "Netral",
            f"{concentration * 1.2:.1f} M" if concentration < 3 else "Cukup",
            "Air atau Etanol disarankan" if solvent not in ['Air', 'Etanol'] else "Cocok"
        ]
    }
    st.table(pd.DataFrame(data))

# Unduh laporan (dummy)
with st.expander("📄 Simpan Laporan"):
    st.download_button(
        label="📥 Unduh Laporan",
        data="Simulasi laporan prediksi reaksi.",
        file_name="laporan_reaksi.txt",
        mime="text/plain"
    )

# Footer
st.caption("🧪 Ini adalah versi ringan tanpa visualisasi molekul.")
