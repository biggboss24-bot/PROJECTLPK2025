import streamlit as st
from PIL import Image
import base64

# --------------------------
# Fungsi untuk latar belakang dari URL
# --------------------------
def add_bg_from_url():
    st.markdown(
        f"""
        <style>
        .stApp {{
            background-image: url("https://images.unsplash.com/photo-1581090700227-1e8e30e106b4");
            background-size: cover;
            background-attachment: fixed;
        }}
        .block-container {{
            background-color: rgba(255, 255, 255, 0.85);
            border-radius: 15px;
            padding: 2rem;
        }}
        </style>
        """,
        unsafe_allow_html=True
    )

# --------------------------
# Tambahkan background
# --------------------------
add_bg_from_url()

# --------------------------
# Judul Aplikasi
# --------------------------
st.markdown("<h1 style='text-align: center; color: navy;'>Reaction Navigator 🔬</h1>", unsafe_allow_html=True)
st.markdown("<h4 style='text-align: center;'>Prediksi Jalur Reaksi Berdasarkan Reaktan yang Diberikan</h4>", unsafe_allow_html=True)

# --------------------------
# Database Reaksi
# --------------------------
reaction_db = {
    "Na + Cl2": "NaCl",
    "H2 + O2": "H2O",
    "C + O2": "CO2",
    "CH4 + 2O2": "CO2 + 2H2O",
    "Fe + S": "FeS",
    "CaCO3": "CaO + CO2",
    "Zn + HCl": "ZnCl2 + H2",
    "HCl + NaOH": "NaCl + H2O",
    "AgNO3 + NaCl": "AgCl + NaNO3",
    "NH3 + HCl": "NH4Cl"
}

# --------------------------
# Input Reaktan
# --------------------------
st.subheader("Masukkan Reaktan")
reactants = st.text_input("Contoh: H2 + O2").strip()

# --------------------------
# Tombol Prediksi
# --------------------------
if st.button("Prediksi Produk"):
    if reactants in reaction_db:
        st.success(f"Produk: **{reaction_db[reactants]}**")
    else:
        st.warning("Reaksi tidak ditemukan dalam database. Silakan coba reaktan lain.")

