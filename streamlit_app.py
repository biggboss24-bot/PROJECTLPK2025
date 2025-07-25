import streamlit as st
import pandas as pd
from datetime import datetime

# 🔢 Inisialisasi memori reaksi selama sesi masih aktif
if 'reaction_history' not in st.session_state:
    st.session_state.reaction_history = []

# 🔹 Judul aplikasi
st.title('Reaction Navigator ')
st.subheader('Prediksi Jalur Reaksi Berdasarkan Struktur Sederhana dan Parameter Reaksi')

# 🛠️ Sidebar - Parameter Reaksi
with st.sidebar:
    st.header('🔧 Parameter Reaksi')
    temperature = st.number_input('Suhu (°C)', min_value=0, max_value=200, value=25, step=1)
    ph = st.slider('pH', 0, 14, 7)
    concentration = st.slider('Konsentrasi (mol/L)', 0.1, 5.0, 1.0)
    solvent = st.selectbox('Pelarut', ['Air', 'Etanol', 'Aseton', 'Dietil Eter', 'DMSO'])
    catalyst_on = st.checkbox('Ada Katalis?')
    catalyst_type = st.selectbox('Jenis Katalis', ['Asam', 'Basa', 'Logam Transisi', 'Enzim']) if catalyst_on else None

# 🧨 Input Reaktan
st.header('🤓 Masukkan Reaktan')
col1, col2 = st.columns(2)
with col1: reactant1 = st.text_input('Reaktan 1 (Contoh: CH3COOH)', 'CH3COOH')
with col2: reactant2 = st.text_input('Reaktan 2 (Contoh: C2H5OH)', 'C2H5OH')

# 🪠 Fungsi Prediksi Jalur Reaksi
def predict_reaction_pathway(r1, r2, temp, ph, conc, solvent, catalyst=None):
    r1 = r1.strip().upper()
    r2 = r2.strip().upper()
    pathways = []

    if ('COOH' in r1 and 'OH' in r2) or ('COOH' in r2 and 'OH' in r1):
        pathways.append({
            'name': 'Esterifikasi Umum',
            'steps': [('Protonasi karbonil', 40), ('Serangan nukleofil', 35), ('Pelepasan air', 25)],
            'products': ['Ester + H2O'],
            'activation_energy': 70 - (ph * 1.5) - (temp * 0.2),
            'thermodynamics': 'Endotermik'
        })
    elif ('HCL' in r1 and 'NAOH' in r2) or ('NAOH' in r1 and 'HCL' in r2) or ('OH' in r1 and 'H' in r2) or ('OH' in r2 and 'H' in r1):
        pathways.append({
            'name': 'Netralisasi Asam-Basa',
            'steps': [('Disosiasi elektrolit', 20), ('Rekombinasi ion', 20), ('Pembentukan air dan garam', 20)],
            'products': ['Garam + H2O'],
            'activation_energy': 45 - (ph * 1.2) - (temp * 0.1),
            'thermodynamics': 'Eksotermik'
        })
    elif 'C=C' in r1 or 'C=C' in r2:
        pathways.append({
            'name': 'Adisi ke Alkena',
            'steps': [('Serangan elektrofilik', 50), ('Intermediate karbokation', 30), ('Adisi nukleofil', 20)],
            'products': ['Produk adisi'],
            'activation_energy': 65 - (ph * 1.2) - (temp * 0.2),
            'thermodynamics': 'Eksotermik'
        })
    elif ('NH2' in r1 and 'COOH' in r2) or ('NH2' in r2 and 'COOH' in r1):
        pathways.append({
            'name': 'Pembentukan Amida',
            'steps': [('Aktivasi asam', 45), ('Serangan NH2', 35), ('Kondensasi', 20)],
            'products': ['Amida + H2O'],
            'activation_energy': 80 - (ph * 1.3) - (temp * 0.15),
            'thermodynamics': 'Endotermik'
        })
    else:
        pathways.append({
            'name': 'Reaksi Umum (Tidak Dikenali)',
            'steps': [('Inisiasi', 40), ('Reaksi berlangsung', 30), ('Produk terbentuk', 20)],
            'products': [f'{r1} + {r2} (prediksi generik)'],
            'activation_energy': 80 - (ph * 1.0) - (temp * 0.15),
            'thermodynamics': 'Tidak Diketahui'
        })

    return pathways

# 🚀 Tombol Prediksi
if st.button('🔍 Prediksi Jalur Reaksi'):
    if not reactant1 or not reactant2:
        st.error('Mohon isi kedua reaktan!')
    else:
        results = predict_reaction_pathway(
            reactant1, reactant2, temperature, ph, concentration, solvent, catalyst_type
        )
        st.success('Prediksi berhasil!')

        for res in results:
            st.subheader(f"🔹 {res['name']}")
            st.markdown(f"**Energi Aktivasi:** {res['activation_energy']:.1f} kJ/mol")
            st.markdown(f"**Termodinamika:** {res['thermodynamics']}")
            st.markdown("### Tahapan Reaksi")
            for i, (step, energy) in enumerate(res['steps'], 1):
                st.write(f"{i}. {step} — {energy} kJ/mol")
            st.markdown("### Produk Prediksi")
            st.code(', '.join(res['products']))
            st.divider()

            st.session_state.reaction_history.append({
                'time': datetime.now().strftime('%H:%M:%S'),
                'r1': reactant1,
                'r2': reactant2,
                'name': res['name'],
                'products': ', '.join(res['products']),
                'Ea': f"{res['activation_energy']:.1f} kJ/mol",
                'thermo': res['thermodynamics']
            })

# 📃 Riwayat Prediksi
st.header("📄 Memori Prediksi Reaksi")
if st.session_state.reaction_history:
    for idx, h in enumerate(reversed(st.session_state.reaction_history), 1):
        st.markdown(f"#### 🔹 Reaksi #{idx} — {h['time']}")
        st.markdown(f"- **Reaktan:** {h['r1']} + {h['r2']}")
        st.markdown(f"- **Jalur:** {h['name']}")
        st.markdown(f"- **Produk:** {h['products']}")
        st.markdown(f"- **Energi Aktivasi:** {h['Ea']}")
        st.markdown(f"- **Termodinamika:** {h['thermo']}")
        st.markdown("---")
else:
    st.info("Belum ada riwayat prediksi.")

# ⚙️ Rekomendasi Kondisi
st.header("⚙️ Rekomendasi Optimasi")
with st.expander("Lihat Saran Optimasi"):
    rekom = {
        'Parameter': ['Suhu', 'pH', 'Konsentrasi', 'Pelarut'],
        'Nilai Sekarang': [f"{temperature} °C", ph, f"{concentration:.1f} M", solvent],
        'Rekomendasi': [
            f"{temperature + 10} °C" if temperature < 100 else "Cukup",
            "pH Asam (3–6)" if ph > 7 else "pH Basa (8–11)" if ph < 6 else "Netral (baik)",
            f"{concentration * 1.2:.1f} M" if concentration < 3 else "Cukup",
            "Air/Etanol disarankan" if solvent not in ['Air', 'Etanol'] else "Sudah cocok"
        ]
    }
    st.table(pd.DataFrame(rekom))

# 📁 Simpan Laporan Sederhana
with st.expander("📁 Simpan Laporan"):
    st.download_button("Unduh Laporan TXT",
        data="Ini adalah hasil prediksi jalur reaksi.",
        file_name="laporan_reaksi.txt",
        mime="text/plain")

st.caption("⚠️ Aplikasi ini adalah simulasi edukatif. Prediksi hanya berdasarkan pola teks, bukan perhitungan struktur molekul.")
