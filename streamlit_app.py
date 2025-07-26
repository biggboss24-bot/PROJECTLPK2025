import streamlit as st
import pandas as pd
from datetime import datetime

# 🌄 Gaya Tampilan (dengan background kimia)
st.markdown("""
    <style>
    .stApp {
        background-color: #001f3f;
        background-image: url("https://images.unsplash.com/photo-1581091012184-5c1d7e0b4a9c?auto=format&fit=crop&w=1200&q=80");
        background-size: cover;
        background-position: center;
        background-attachment: fixed;
        color: white;
    }
    .block-container {
        background-color: rgba(0, 0, 0, 0.6);
        padding: 2rem;
        border-radius: 12px;
    }
    h1, h2, h3 {
        color: cyan;
    }
    </style>
""", unsafe_allow_html=True)

# 📚 Database Reaksi Kimia
REACTION_DATABASE = [
    {
        'reactants': ['NaOH', 'HCl'],
        'name': 'Netralisasi',
        'products': 'NaCl + H₂O',
        'Ea': '20 kJ/mol',
        'thermo': 'Eksotermik'
    },
    {
        'reactants': ['CH₃COOH', 'C₂H₅OH'],
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
        'reactants': ['CaCO₃', 'HCl'],
        'name': 'Reaksi Asam-Karbonat',
        'products': 'CaCl₂ + CO₂ + H₂O',
        'Ea': '35 kJ/mol',
        'thermo': 'Eksotermik'
    }
]

# 🔍 Fungsi Pencarian Reaksi
def cari_reaksi(r1, r2):
    r1 = r1.strip().lower()
    r2 = r2.strip().lower()
    for item in REACTION_DATABASE:
        reaktan_db = [r.lower() for r in item['reactants']]
        if sorted([r1, r2]) == sorted(reaktan_db):
            return item
    return None

# 🔁 Inisialisasi Riwayat
if 'history' not in st.session_state:
    st.session_state.history = []

# 🧪 Judul
st.title("🔬 Reaction Navigator")
st.subheader("Prediksi Jalur Reaksi Berdasarkan Dua Reaktan")

# 📥 Input Reaktan
st.markdown("## Masukkan Reaktan")
reaktan1 = st.text_input("Reaktan 1")
reaktan2 = st.text_input("Reaktan 2")

# 🔘 Tombol Prediksi
if st.button("Prediksi Reaksi"):
    if reaktan1 and reaktan2:
        hasil = cari_reaksi(reaktan1, reaktan2)
        if hasil:
            st.success("✅ Reaksi ditemukan!")
            st.write(f"**Nama Reaksi:** {hasil['name']}")
            st.write(f"**Produk:** {hasil['products']}")
            st.write(f"**Energi Aktivasi:** {hasil['Ea']}")
            st.write(f"**Sifat Termodinamika:** {hasil['thermo']}")
            
            st.session_state.history.append({
                'waktu': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                'r1': reaktan1,
                'r2': reaktan2,
                'hasil': hasil
            })
        else:
            st.warning("❌ Reaksi tidak ditemukan dalam database.")
    else:
        st.warning("⚠️ Masukkan kedua reaktan terlebih dahulu.")

# 📜 Tampilkan Riwayat Reaksi
if st.session_state.history:
    st.markdown("## 📚 Riwayat Prediksi Reaksi")
    for i, h in enumerate(reversed(st.session_state.history), 1):
        st.markdown(f"### Reaksi #{i}")
        st.markdown(f"- **Waktu:** {h['waktu']}")
        st.markdown(f"- **Reaktan:** {h['r1']} + {h['r2']}")
        st.markdown(f"- **Produk:** {h['hasil']['products']}")
        st.markdown(f"- **Nama Reaksi:** {h['hasil']['name']}")
        st.markdown(f"- **Ea:** {h['hasil']['Ea']}")
        st.markdown(f"- **Termodinamika:** {h['hasil']['thermo']}")
        st.markdown("---")

# 💾 Unduh Laporan
with st.expander("📥 Unduh Laporan Prediksi"):
    if st.session_state.history:
        isi_laporan = ""
        for i, h in enumerate(st.session_state.history, 1):
            isi_laporan += f"Reaksi #{i} ({h['waktu']})\n"
            isi_laporan += f"Reaktan: {h['r1']} + {h['r2']}\n"
            isi_laporan += f"Produk: {h['hasil']['products']}\n"
            isi_laporan += f"Nama Reaksi: {h['hasil']['name']}\n"
            isi_laporan += f"Energi Aktivasi: {h['hasil']['Ea']}\n"
            isi_laporan += f"Sifat Termodinamika: {h['hasil']['thermo']}\n"
            isi_laporan += "-"*40 + "\n"
        st.download_button("📄 Unduh Laporan", isi_laporan, "laporan_reaksi.txt", "text/plain")
    else:
        st.info("Belum ada reaksi yang diprediksi.")
