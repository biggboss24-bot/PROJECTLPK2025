import streamlit as st
import pandas as pd

st.title('🧪 Reaction Pathway Predictor (Versi Ringan)')
st.subheader('Prediksi Jalur Reaksi Berdasarkan Parameter dan Reaktan')

# Sidebar input parameter
with st.sidebar:
    st.header('🔧 Parameter Reaksi')
    temperature = st.slider('Suhu (°C)', 0, 200, 25)
    ph = st.slider('pH', 0, 14, 7)
    concentration = st.slider('Konsentrasi (mol/L)', 0.1, 5.0, 1.0)
    solvent = st.selectbox('Pelarut', ['Air', 'Etanol', 'Aseton', 'Dietil Eter', 'DMSO'])
    catalyst = st.checkbox('Ada Katalis?')
    if catalyst:
        catalyst_type = st.selectbox('Jenis Katalis', ['Asam', 'Basa', 'Logam Transisi', 'Enzim'])

# Input reaktan
st.header('🔬 Masukkan Reaktan')
col1, col2 = st.columns(2)
with col1:
    reactant1 = st.text_input('Reaktan 1 (contoh: HCl)', 'C=O')
with col2:
    reactant2 = st.text_input('Reaktan 2 (contoh: NaOH)', 'O')

# Fungsi prediksi
def predict_reaction_pathway(r1, r2, temp, ph, conc, solvent, catalyst=None):
    r1 = r1.strip()
    r2 = r2.strip()

    pathways = []

    if (r1 == 'C=O' and r2 == 'O') or (r1 == 'O' and r2 == 'C=O'):
        pathways.append({
            'name': 'Hidrasi Aldehid',
            'steps': [
                ('Pembentukan intermediate tetrahedral', 50),
                ('Proton transfer', 30),
                ('Pembentukan hidrat', 20)
            ],
            'products': ['OC(O)C'],
            'activation_energy': 75 - (ph * 2) - (temp * 0.1),
            'thermodynamics': 'Eksotermik'
        })
    elif (r1 == 'CC(=O)O' and r2 == 'CO') or (r1 == 'CO' and r2 == 'CC(=O)O'):
        pathways.append({
            'name': 'Esterifikasi',
            'steps': [
                ('Protonasi karbonil', 60),
                ('Serangan nukleofilik', 40),
                ('Dehidrasi', 30)
            ],
            'products': ['CC(=O)OC + H₂O'],
            'activation_energy': 90 - (ph * 3) - (temp * 0.2),
            'thermodynamics': 'Endotermik'
        })
    elif (r1 == 'HCl' and r2 == 'NaOH') or (r1 == 'NaOH' and r2 == 'HCl'):
        pathways.append({
            'name': 'Netralisasi Asam-Basa',
            'steps': [
                ('Disosiasi HCl dan NaOH', 20),
                ('Pembentukan ion Na⁺ dan OH⁻', 25),
                ('Pembentukan NaCl dan H₂O', 15)
            ],
            'products': ['NaCl + H₂O'],
            'activation_energy': 45 - (ph * 1.2) - (temp * 0.1),
            'thermodynamics': 'Eksotermik'
        })
    else:
        pathways.append({
            'name': 'Reaksi Adisi Umum',
            'steps': [
                ('Inisiasi', 40),
                ('Propagasi', 30),
                ('Terminasi', 20)
            ],
            'products': [f'{r1} + {r2} (produk tidak dikenali)'],
            'activation_energy': 85 - (ph * 1.5) - (temp * 0.15),
            'thermodynamics': 'Eksotermik'
        })

    return pathways

# Tombol prediksi
if st.button('🔍 Prediksi Jalur Reaksi'):
    if not reactant1 or not reactant2:
        st.error('❗ Mohon masukkan kedua reaktan!')
    else:
        with st.spinner('⏳ Menganalisis jalur reaksi...'):
            results = predict_reaction_pathway(
                reactant1, reactant2, temperature, ph, concentration, solvent, catalyst
            )

            st.success('✅ Prediksi berhasil!')

            for path in results:
                st.subheader(f"🔷 Jalur Reaksi: {path['name']}")
                st.markdown(f"**Energi Aktivasi:** {path['activation_energy']:.1f} kJ/mol")
                st.markdown(f"**Termodinamika:** {path['thermodynamics']}")

                st.markdown("### 🧪 Tahapan Reaksi")
                for i, (step, energy) in enumerate(path['steps'], 1):
                    st.write(f"{i}. {step} — {energy} kJ/mol")

                st.markdown("### 🧬 Produk Prediksi")
                for prod in path['products']:
                    st.code(prod)

                st.divider()

# Rekomendasi kondisi reaksi
st.header("⚙️ Rekomendasi Optimasi")
with st.expander("🔁 Saran Kondisi yang Lebih Efektif"):
    rekomendasi = {
        'Parameter': ['Suhu', 'pH', 'Konsentrasi', 'Pelarut'],
        'Nilai Saat Ini': [f"{temperature}°C", ph, f"{concentration:.1f} M", solvent],
        'Rekomendasi': [
            f"{temperature + 10}°C" if temperature < 100 else "Cukup",
            "Asam (pH 3–6)" if ph > 7 else "Basa (pH 8–11)" if ph < 6 else "Netral (baik)",
            f"{concentration * 1.2:.1f} M" if concentration < 3 else "Cukup",
            "Air atau Etanol disarankan" if solvent not in ['Air', 'Etanol'] else "Sudah sesuai"
        ]
    }
    st.table(pd.DataFrame(rekomendasi))

# Simpan hasil analisis (dummy)
with st.expander("📁 Simpan Laporan"):
    st.download_button(
        label="📥 Unduh Laporan TXT",
        data="Ini contoh hasil simulasi reaksi kimia. Laporan lengkap bisa dikembangkan.",
        file_name="laporan_reaksi.txt",
        mime="text/plain"
    )

st.caption("⚠️ Aplikasi ini adalah simulasi edukatif. Hasil eksperimen nyata bisa berbeda.")
