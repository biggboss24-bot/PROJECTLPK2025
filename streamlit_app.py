import streamlit as st
import pandas as pd
from datetime import datetime

# ───────────────────────────────────────────────────────────
# ⚙️  Inisialisasi memori reaksi selama sesi masih hidup
# ───────────────────────────────────────────────────────────
if 'reaction_history' not in st.session_state:
    st.session_state.reaction_history = []

# ───────────────────────────────────────────────────────────
# 🏷️  Judul & sub-header
# ───────────────────────────────────────────────────────────
st.title('🧪 Reaction Pathway Predictor (Versi Ringan + Memori)')
st.subheader('Prediksi Jalur Reaksi Berdasarkan Parameter dan Reaktan')

# ───────────────────────────────────────────────────────────
# 🛠️  Sidebar – parameter reaksi
# ───────────────────────────────────────────────────────────
with st.sidebar:
    st.header('🔧 Parameter Reaksi')
    temperature = st.slider('Suhu (°C)', 0, 200, 25)
    ph          = st.slider('pH', 0, 14, 7)
    concentration = st.slider('Konsentrasi (mol/L)', 0.1, 5.0, 1.0)
    solvent     = st.selectbox('Pelarut', ['Air', 'Etanol', 'Aseton', 'Dietil Eter', 'DMSO'])
    catalyst_on = st.checkbox('Ada Katalis?')
    catalyst_type = None
    if catalyst_on:
        catalyst_type = st.selectbox('Jenis Katalis', ['Asam', 'Basa', 'Logam Transisi', 'Enzim'])

# ───────────────────────────────────────────────────────────
# 🧬 Input reaktan
# ───────────────────────────────────────────────────────────
st.header('🔬 Masukkan Reaktan')
c1, c2 = st.columns(2)
with c1: reactant1 = st.text_input('Reaktan 1', 'C=O')
with c2: reactant2 = st.text_input('Reaktan 2', 'O')

# ───────────────────────────────────────────────────────────
# 🔮 Fungsi prediksi jalur reaksi
# ───────────────────────────────────────────────────────────
def predict_reaction_pathway(r1, r2, temp, ph, conc, solvent, catalyst=None):
    r1, r2 = r1.strip(), r2.strip()
    paths = []

    # Hidrasi aldehid
    if {r1, r2} == {'C=O', 'O'}:
        paths.append(dict(
            name='Hidrasi Aldehid',
            steps=[('Intermediate tetrahedral', 50),
                   ('Proton transfer', 30),
                   ('Hidrat terbentuk', 20)],
            products=['OC(O)C'],
            activation_energy=75 - (ph*2) - (temp*0.1),
            thermodynamics='Eksotermik'
        ))

    # Esterifikasi sederhana
    elif {r1, r2} == {'CC(=O)O', 'CO'}:
        paths.append(dict(
            name='Esterifikasi',
            steps=[('Protonasi karbonil', 60),
                   ('Serangan nukleofilik', 40),
                   ('Dehidrasi', 30)],
            products=['CC(=O)OC + H₂O'],
            activation_energy=90 - (ph*3) - (temp*0.2),
            thermodynamics='Endotermik'
        ))

    # Netralisasi HCl + NaOH
    elif {r1, r2} == {'HCl', 'NaOH'}:
        paths.append(dict(
            name='Netralisasi Asam-Basa',
            steps=[('Disosiasi elektrolit', 20),
                   ('Pembentukan ion Na⁺ & OH⁻', 25),
                   ('Pembentukan NaCl & H₂O', 15)],
            products=['NaCl + H₂O'],
            activation_energy=45 - (ph*1.2) - (temp*0.1),
            thermodynamics='Eksotermik'
        ))

    # Fallback
    else:
        paths.append(dict(
            name='Reaksi Adisi Umum',
            steps=[('Inisiasi', 40),
                   ('Propagasi', 30),
                   ('Terminasi', 20)],
            products=[f'{r1} + {r2} (produk tidak dikenali)'],
            activation_energy=85 - (ph*1.5) - (temp*0.15),
            thermodynamics='Eksotermik'
        ))

    return paths

# ───────────────────────────────────────────────────────────
# 🚀 Tombol prediksi
# ───────────────────────────────────────────────────────────
if st.button('🔍 Prediksi Jalur Reaksi'):
    if not reactant1 or not reactant2:
        st.error('❗ Mohon isi kedua reaktan.')
    else:
        with st.spinner('Menganalisis...'):
            results = predict_reaction_pathway(
                reactant1, reactant2, temperature, ph, concentration, solvent, catalyst_type
            )

        st.success('✅ Prediksi berhasil!')

        # ── Tampilkan & simpan ke memori
        for res in results:
            st.subheader(f"🔷 {res['name']}")
            st.markdown(f"**Energi Aktivasi:** {res['activation_energy']:.1f} kJ/mol")
            st.markdown(f"**Termodinamika:** {res['thermodynamics']}")

            st.markdown("### 🧪 Tahapan Reaksi")
            for i, (step, en) in enumerate(res['steps'], 1):
                st.write(f"{i}. {step} — {en} kJ/mol")

            st.markdown("### 🧬 Produk Prediksi")
            st.code(', '.join(res['products']))

            st.divider()

            # Simpan tiap hasil ke session_state
            st.session_state.reaction_history.append(dict(
                time=datetime.now().strftime('%H:%M:%S'),
                r1=reactant1,
                r2=reactant2,
                name=res['name'],
                products=', '.join(res['products']),
                Ea=f"{res['activation_energy']:.1f} kJ/mol",
                thermo=res['thermodynamics']
            ))

# ───────────────────────────────────────────────────────────
# 📚 Memori Prediksi Reaksi
# ───────────────────────────────────────────────────────────
st.header("🧾 Memori Prediksi Reaksi (sesi ini)")
if st.session_state.reaction_history:
    for idx, h in enumerate(reversed(st.session_state.reaction_history), 1):
        st.markdown(f"#### 🔹 Reaksi #{idx} — {h['time']}")
        st.markdown(f"* **Reaktan:** {h['r1']} + {h['r2']}")
        st.markdown(f"* **Jalur:** {h['name']}")
        st.markdown(f"* **Produk:** {h['products']}")
        st.markdown(f"* **Energi Aktivasi:** {h['Ea']}")
        st.markdown(f"* **Termodinamika:** {h['thermo']}")
        st.markdown("---")
else:
    st.info("Belum ada riwayat prediksi.")

# ───────────────────────────────────────────────────────────
# ⚙️ Rekomendasi optimasi
# ───────────────────────────────────────────────────────────
st.header("⚙️ Rekomendasi Optimasi")
with st.expander("Lihat saran kondisi yang lebih efektif"):
    rec = {
        'Parameter'      : ['Suhu', 'pH', 'Konsentrasi', 'Pelarut'],
        'Nilai Saat Ini' : [f"{temperature} °C", ph, f"{concentration:.1f} M", solvent],
        'Rekomendasi'    : [
            f"{temperature+10} °C" if temperature < 100 else "Cukup",
            "Asam (3-6)" if ph > 7 else "Basa (8-11)" if ph < 6 else "Netral",
            f"{concentration*1.2:.1f} M" if concentration < 3 else "Cukup",
            "Air/Etanol disarankan" if solvent not in ['Air', 'Etanol'] else "Sudah sesuai"
        ]
    }
    st.table(pd.DataFrame(rec))

# ───────────────────────────────────────────────────────────
# 💾 Unduh laporan dummy
# ───────────────────────────────────────────────────────────
with st.expander("📥 Simpan Laporan Sederhana (.txt)"):
    st.download_button(
        "Unduh",
        data="Ini contoh laporan simulasi reaksi kimia.",
        file_name="laporan_reaksi.txt",
        mime="text/plain"
    )

st.caption("ℹ️ Aplikasi simulasi edukatif—hasil eksperimen nyata dapat berbeda.")
