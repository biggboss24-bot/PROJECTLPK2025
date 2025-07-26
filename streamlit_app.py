import streamlit as st
import pandas as pd
from datetime import datetime

# 🌄 Gaya tampilan latar belakang
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
    h1, h2, h3 {
        color: #00ffff;
    }
    </style>
""", unsafe_allow_html=True)

# 🌟 Database awal reaksi
if "reaction_database" not in st.session_state:
    st.session_state.reaction_database = [
        {
            'reactants': ['NaOH', 'HCl'],
            'name': 'Netralisasi',
            'products': 'NaCl + H₂O',
            'Ea': '20 kJ/mol',
            'thermo': 'Eksotermik',
            'type': 'Netralisasi'
        },
        {
            'reactants': ['CH₄', 'O₂'],
            'name': 'Pembakaran Metana',
            'products': 'CO₂ + H₂O',
            'Ea': '60 kJ/mol',
            'thermo': 'Eksotermik',
            'type': 'Pembakaran'
        }
    ]

# 🕘 Inisialisasi riwayat reaksi
if 'reaction_history' not in st.session_state:
    st.session_state.reaction_history = []

# 🔎 Fungsi pencarian reaksi
def cari_reaksi(reaktan1, reaktan2):
    r1 = reaktan1.strip().lower()
    r2 = reaktan2.strip().lower()
    for reaksi in st.session_state.reaction_database:
        reactant_set = set([r.lower() for r in reaksi['reactants']])
        if set([r1, r2]) == reactant_set:
            return reaksi
    return None

# 🧪 Judul utama
st.title('🔬 Reaction Navigator')
st.subheader('Prediksi Jalur Reaksi Kimia & Penambahan Reaksi')

# 👉 Input reaktan untuk prediksi
st.markdown("## 🔍 Prediksi Reaksi")
reaktan1 = st.text_input('Reaktan 1')
reaktan2 = st.text_input('Reaktan 2')

if st.button("Prediksi Jalur Reaksi"):
    hasil = cari_reaksi(reaktan1, reaktan2)
    if hasil:
        st.success(f"Reaksi Ditemukan: {hasil['name']}")
        st.write(f"**Produk:** {hasil['products']}")
        st.write(f"**Jenis:** {hasil['type']}")
        st.write(f"**Energi Aktivasi:** {hasil['Ea']}")
        st.write(f"**Termodinamika:** {hasil['thermo']}")
        # Simpan ke riwayat
        st.session_state.reaction_history.append({
            'time': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'r1': reaktan1,
            'r2': reaktan2,
            'name': hasil['name'],
            'products': hasil['products'],
            'type': hasil['type'],
            'Ea': hasil['Ea'],
            'thermo': hasil['thermo']
        })
    else:
        st.warning("Reaksi tidak ditemukan di database!")

# ➕ Tambah reaksi ke database
st.markdown("## ➕ Tambah Reaksi Baru")
with st.form("tambah_reaksi"):
    new_r1 = st.text_input("Reaktan 1", key="new_r1")
    new_r2 = st.text_input("Reaktan 2", key="new_r2")
    new_product = st.text_input("Produk Reaksi")
    new_type = st.selectbox("Jenis Reaksi", ['Netralisasi', 'Pembakaran', 'Dekomposisi', 'Sintesis', 'Substitusi'])
    new_name = st.text_input("Nama Reaksi")
    new_Ea = st.text_input("Energi Aktivasi (misal: 35 kJ/mol)")
    new_thermo = st.selectbox("Sifat Termodinamika", ['Eksotermik', 'Endotermik', 'Netral'])

    submitted = st.form_submit_button("✅ Tambahkan ke Database")
    if submitted:
        st.session_state.reaction_database.append({
            'reactants': [new_r1, new_r2],
            'name': new_name,
            'products': new_product,
            'Ea': new_Ea,
            'thermo': new_thermo,
            'type': new_type
        })
        st.success("Reaksi berhasil ditambahkan!")

# 🧾 Tampilkan riwayat
if st.session_state.reaction_history:
    st.markdown("## 📜 Riwayat Prediksi")
    for i, h in enumerate(reversed(st.session_state.reaction_history), 1):
        st.markdown(f"### 🔁 Reaksi #{i}")
        st.markdown(f"- ⏰ Waktu: {h['time']}")
        st.markdown(f"- 🧪 Reaktan: {h['r1']} + {h['r2']}")
        st.markdown(f"- 🔬 Nama Reaksi: {h['name']}")
        st.markdown(f"- ⚗️ Produk: {h['products']}")
        st.markdown(f"- 🧭 Jenis: {h['type']}")
        st.markdown(f"- 🔥 Energi Aktivasi: {h['Ea']}")
        st.markdown(f"- 🌡️ Termodinamika: {h['thermo']}")
        st.markdown("---")
