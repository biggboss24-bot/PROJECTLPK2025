import streamlit as st
import pandas as pd
from datetime import datetime
import base64

# --- Styling: Gambar Background + Font ---
def set_background(image_file):
    with open(image_file, "rb") as f:
        data = base64.b64encode(f.read()).decode()
    st.markdown(f"""
        <style>
        .stApp {{
            background-image: url("data:image/jpg;base64,{data}");
            background-size: cover;
            background-attachment: fixed;
        }}
        .block-container {{
            background-color: rgba(0,0,0,0.5);
            padding: 2rem;
            border-radius: 10px;
        }}
        h1, h2, h3, p, span {{
            color: white;
        }}
        </style>
    """, unsafe_allow_html=True)

set_background("lab_bg.jpg")

# --- Inisialisasi session_state untuk database dan riwayat ---
if 'reaction_history' not in st.session_state:
    st.session_state.reaction_history = []

if 'custom_reactions' not in st.session_state:
    st.session_state.custom_reactions = {
        'Netralisasi': {'Asam + Basa': 'Garam + Air'},
        'Pembakaran': {'CH4 + O2': 'CO2 + H2O'},
        'Dekomposisi': {'H2O2': 'H2O + O2'}
    }

# --- Judul Aplikasi ---
st.title("🔬 Reaction Navigator")
st.subheader("Prediksi Jalur Reaksi Berdasarkan Parameter dan Jenis Reaksi")

# --- Sidebar: Parameter Reaksi ---
with st.sidebar:
    st.header("🔧 Parameter Reaksi")
    jenis = st.selectbox("Jenis Reaksi", ["Netralisasi", "Pembakaran", "Dekomposisi"])
    suhu = st.slider("Suhu (°C)", 0, 300, 25)
    ph = st.slider("pH", 0, 14, 7)
    konsentrasi = st.number_input("Konsentrasi (mol/L)", 0.1, 5.0, 1.0)

# --- Input Reaktan ---
st.markdown("### 🧪 Input Reaktan")
reaktan1 = st.text_input("Reaktan 1")
reaktan2 = st.text_input("Reaktan 2 (jika ada)")

# --- Prediksi Reaksi ---
if st.button("🔍 Prediksi Reaksi"):
    key = f"{reaktan1} + {reaktan2}" if reaktan2 else reaktan1
    hasil = st.session_state.custom_reactions.get(jenis, {}).get(key, None)

    if hasil:
        st.success("Prediksi reaksi ditemukan!")
        st.markdown(f"**Reaksi:** {key} → {hasil}")
        st.image("https://upload.wikimedia.org/wikipedia/commons/thumb/1/13/Combustion_Reaction.svg/800px-Combustion_Reaction.svg.png", width=400)

        st.session_state.reaction_history.append({
            'waktu': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'jenis': jenis,
            'reaksi': f"{key} → {hasil}"
        })
    else:
        st.error("Reaksi belum ada di database. Tambahkan di bawah!")

# --- Tambah Reaksi ke Database ---
st.markdown("### ➕ Tambahkan Reaksi ke Database")
with st.form("form_tambah"):
    jenis_tambah = st.selectbox("Jenis Reaksi Baru", ["Netralisasi", "Pembakaran", "Dekomposisi"])
    reaktan_baru = st.text_input("Reaktan (cth: HCl + NaOH)")
    produk_baru = st.text_input("Produk (cth: NaCl + H2O)")
    tambah = st.form_submit_button("Tambah Reaksi")
    if tambah and reaktan_baru and produk_baru:
        if jenis_tambah not in st.session_state.custom_reactions:
            st.session_state.custom_reactions[jenis_tambah] = {}
        st.session_state.custom_reactions[jenis_tambah][reaktan_baru] = produk_baru
        st.success("Reaksi berhasil ditambahkan!")

# --- Riwayat Reaksi ---
if st.session_state.reaction_history:
    st.markdown("### 📜 Riwayat Reaksi")
    for idx, h in enumerate(st.session_state.reaction_history[::-1], 1):
        st.markdown(f"**{idx}. {h['waktu']}** — *{h['jenis']}*  \n🧪 {h['reaksi']}")

# --- Unduh Riwayat ---
if st.session_state.reaction_history:
    if st.download_button("⬇️ Unduh Riwayat Reaksi", 
                          data="\n".join([f"{h['waktu']} - {h['reaksi']}" for h in st.session_state.reaction_history]), 
                          file_name="riwayat_reaksi.txt"):
        st.success("Berhasil diunduh!")
