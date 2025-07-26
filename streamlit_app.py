import streamlit as st
import pandas as pd
from datetime import datetime

# 🌄 Gaya tampilan
st.markdown("""
    <style>
    .stApp {
        background-color: #001f3f;
        background-image: url("https://images.unsplash.com/photo-1581091012184-5c1d7e0b4a9c?auto=format&fit=crop&w=1200&q=80");
        background-size: cover;
        background-position: center;
        background-attachment: fixed;
        color: #ffffff;
    }
    .block-container {
        background-color: rgba(0, 0, 0, 0.6);
        padding: 2rem;
        border-radius: 10px;
    }
    h1, h2, h3, h4 {
        color: #00ffff;
    }
    </style>
""", unsafe_allow_html=True)

# 🧪 Database reaksi kimia
REACTION_DATABASE = [
    {
        'reactants': ['NaOH', 'HCl'],
        'name': 'Reaksi Netralisasi',
        'products': 'NaCl + H₂O',
        'Ea': '20 kJ/mol',
        'thermo': 'Eksotermik'
    },
    {
        'reactants': ['C₂H₅OH', 'CH₃COOH'],
        'name': 'Esterifikasi',
        'products': 'CH₃COOC₂H₅ + H₂O',
        'Ea': '65 kJ/mol',
        'thermo': 'Endotermik'
    },
    {
        'reactants': ['AgNO₃', 'NaCl'],
        'name': 'Presipitasi',
        'products': 'AgCl↓ + NaNO₃',
        'Ea': '10 kJ/mol',
        'thermo': 'Netral'
    },
    {
        'reactants': ['H₂', 'Cl₂'],
        'name': 'Reaksi Sintesis',
        'products': '2HCl',
        'Ea': '75 kJ/mol',
        'thermo': 'Eksotermik'
    },
    {
        'reactants': ['CaCO₃', 'HCl'],
        'name': 'Reaksi Asam Basa',
        'products': 'CaCl₂ + CO₂ + H₂O',
        'Ea': '35 kJ/mol',
        'thermo': 'Eksotermik'
    }
]

# 🔎 Fungsi pencarian reaksi
def cari_reaksi(reaktan1, reaktan2):
    r1 = reaktan1.strip().lower()
    r2 = reaktan2.strip().lower()
    for reaksi in REACTION_DATABASE:
        reactant_set = set([r.lower() for r in reaksi['reactants']])
        if set([r1, r2]) == reactant_set:
            return reaksi
    return None

# 🔢 Inisialisasi session
if 'reaction_history' not in st.session_state:
    st.session_state.reaction_history = []

# 🧪 Judul
st.title('🔬 Reaction Navigator')
st.subheader('Prediksi Jalur Reaksi Berdasarkan Input Senyawa')

# 🗞 Sidebar parameter
with st.sidebar:
    st.header('Parameter Reaksi')
    temperature = st.number_input('Suhu (°C)', 0, 200, 25)
    ph = st.slider('pH', 0, 14, 7)
    konsentrasi = st.slider('Konsentrasi (mol/L)', 0.1, 5.0, 1.0)
    pelarut = st.selectbox('Pelarut', ['Air', 'Etanol', 'Aseton', 'DMSO'])
    katalis = st.checkbox('Menggunakan Katalis?')
    katalis_jenis = st.selectbox('Jenis Katalis', ['Asam', 'Basa', 'Enzim']) if katalis else None

# 🗞 Input reaktan
st.markdown('## Masukkan Reaktan')
reaktan1 = st.text_input('Reaktan 1')
reaktan2 = st.text_input('Reaktan 2')

# 🚀 Tombol prediksi
if st.button('Prediksi Jalur Reaksi'):
    if reaktan1 and reaktan2:
        hasil = cari_reaksi(reaktan1, reaktan2)
        if hasil:
            st.success('Jalur reaksi ditemukan!')
            st.write(f"**Nama Reaksi:** {hasil['name']}")
            st.write(f"**Produk:** {hasil['products']}")
            st.write(f"**Energi Aktivasi:** {hasil['Ea']}")
            st.write(f"**Termodinamika:** {hasil['thermo']}")

            # Simpan ke riwayat
            st.session_state.reaction_history.append({
                'time': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                'r1': reaktan1,
                'r2': reaktan2,
                'name': hasil['name'],
                'products': hasil['products'],
                'Ea': hasil['Ea'],
                'thermo': hasil['thermo']
            })
        else:
            st.warning('Reaksi tidak ditemukan dalam database.')
    else:
        st.warning('Harap masukkan kedua reaktan.')

# 🗞 Riwayat prediksi
if st.session_state.reaction_history:
    st.markdown('## Riwayat Prediksi')
    for i, h in enumerate(reversed(st.session_state.reaction_history), 1):
        st.markdown(f"### Reaksi #{i}")
        st.markdown(f"- Waktu: {h['time']}")
        st.markdown(f"- Reaktan: {h['r1']} + {h['r2']}")
        st.markdown(f"- Jalur: {h['name']}")
        st.markdown(f"- Produk: {h['products']}")
        st.markdown(f"- Energi Aktivasi: {h['Ea']}")
        st.markdown(f"- Termodinamika: {h['thermo']}")
        st.markdown("---")

# 📂 Ekspor laporan
with st.expander("📄 Unduh Laporan Reaksi"):
    if st.session_state.reaction_history:
        laporan = ""
        for i, h in enumerate(st.session_state.reaction_history, 1):
            laporan += f"Reaksi #{i} - {h['time']}\n"
            laporan += f"Reaktan: {h['r1']} + {h['r2']}\n"
            laporan += f"Jalur: {h['name']}\n"
            laporan += f"Produk: {h['products']}\n"
            laporan += f"Energi Aktivasi: {h['Ea']}\n"
            laporan += f"Termodinamika: {h['thermo']}\n"
            laporan += "-"*40 + "\n"
        st.download_button("📅 Unduh TXT",
                           data=laporan,
                           file_name="laporan_reaksi.txt",
                           mime="text/plain")
    else:
        st.info("Belum ada reaksi untuk disimpan.")
