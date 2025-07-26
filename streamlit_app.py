import streamlit as st
import pandas as pd
from datetime import datetime
import base64

# 🌄 Tambahkan Background dengan Gambar Laboratorium
def set_background():
    st.markdown(
        f"""
        <style>
        .stApp {{
            background-image: url('https://i.imgur.com/2nCt3Sbl.jpg');
            background-size: cover;
            background-position: center;
            background-attachment: fixed;
            color: #ffffff;
        }}
        .block-container {{
            background-color: rgba(0, 0, 0, 0.5);
            padding: 2rem;
            border-radius: 10px;
        }}
        h1, h2, h3, h4 {{
            color: #00e5ff;
            text-shadow: 1px 1px 4px #000000;
        }}
        p, li, span, div {{
            color: #f8f9fa;
        }}
        </style>
        """,
        unsafe_allow_html=True
    )

set_background()

# 🔢 Inisialisasi memori reaksi selama sesi masih aktif
if 'reaction_history' not in st.session_state:
    st.session_state.reaction_history = []

# Contoh database reaksi sederhana
reaction_db = {
    ('HCl', 'NaOH'): {
        'name': 'Netralisasi Asam-Basa',
        'products': 'NaCl + H₂O',
        'Ea': '10 kJ/mol',
        'thermo': 'Eksotermik'
    },
    ('CH3Br', 'NaOH'): {
        'name': 'Substitusi Nukleofilik SN2',
        'products': 'CH3OH + NaBr',
        'Ea': '50 kJ/mol',
        'thermo': 'Eksotermik'
    },
    ('CH2=CH2', 'HBr'): {
        'name': 'Adisi Elektrofilik',
        'products': 'CH3CH2Br',
        'Ea': '60 kJ/mol',
        'thermo': 'Eksotermik'
    },
}

# Judul aplikasi
st.title('Reaction Navigator')
st.subheader('Prediksi Jalur Reaksi Berdasarkan Struktur Sederhana dan Parameter Reaksi')

# Sidebar - Parameter Reaksi
with st.sidebar:
    st.header('Parameter Reaksi')
    temperature = st.number_input('Suhu (°C)', min_value=0, max_value=200, value=25, step=1)
    ph = st.slider('pH', 0, 14, 7)
    concentration = st.slider('Konsentrasi (mol/L)', 0.1, 5.0, 1.0)
    solvent = st.selectbox('Pelarut', ['Air', 'Etanol', 'Aseton', 'Dietil Eter', 'DMSO'])
    catalyst_on = st.checkbox('Ada Katalis?')
    catalyst_type = st.selectbox('Jenis Katalis', ['Asam', 'Basa', 'Logam Transisi', 'Enzim']) if catalyst_on else None

# Input senyawa
st.markdown('### Input Senyawa')
reactant_1 = st.text_input('Reaktan 1', '')
reactant_2 = st.text_input('Reaktan 2', '')

# Prediksi jalur reaksi
if st.button('Prediksi Jalur Reaksi'):
    if reactant_1 and reactant_2:
        key = (reactant_1.strip(), reactant_2.strip())
        reverse_key = (reactant_2.strip(), reactant_1.strip())

        if key in reaction_db:
            predicted_pathway = reaction_db[key]
        elif reverse_key in reaction_db:
            predicted_pathway = reaction_db[reverse_key]
        else:
            predicted_pathway = {
                'name': 'Reaksi Tidak Dikenali',
                'products': f'{reactant_1}-{reactant_2}',
                'Ea': 'N/A',
                'thermo': 'Tidak diketahui'
            }

        st.success('Jalur reaksi berhasil diprediksi!')
        st.write(f"**Nama Jalur Reaksi:** {predicted_pathway['name']}")
        st.write(f"**Produk:** {predicted_pathway['products']}")
        st.write(f"**Energi Aktivasi:** {predicted_pathway['Ea']}")
        st.write(f"**Termodinamika:** {predicted_pathway['thermo']}")

        # Simpan riwayat
        st.session_state.reaction_history.append({
            'time': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'r1': reactant_1,
            'r2': reactant_2,
            'name': predicted_pathway['name'],
            'products': predicted_pathway['products'],
            'Ea': predicted_pathway['Ea'],
            'thermo': predicted_pathway['thermo']
        })
    else:
        st.warning('Mohon masukkan kedua reaktan.')

# Tampilkan riwayat
if st.session_state.reaction_history:
    st.markdown('### Riwayat Prediksi')
    for idx, h in enumerate(st.session_state.reaction_history[::-1], 1):
        st.markdown(f"**Reaksi #{len(st.session_state.reaction_history) - idx + 1}:**")
        st.markdown(f"- Waktu: {h['time']}")
        st.markdown(f"- Reaktan: {h['r1']} + {h['r2']}")
        st.markdown(f"- Jalur: {h['name']}")
        st.markdown(f"- Produk: {h['products']}")
        st.markdown(f"- Energi Aktivasi: {h['Ea']}")
        st.markdown(f"- Termodinamika: {h['thermo']}")
        st.markdown('---')

# Simpan laporan
with st.expander("Simpan Laporan"):
    if st.session_state.reaction_history:
        laporan = ""
        for idx, h in enumerate(st.session_state.reaction_history, 1):
            laporan += f"Reaksi #{idx} — {h['time']}\n"
            laporan += f"Reaktan: {h['r1']} + {h['r2']}\n"
            laporan += f"Jalur: {h['name']}\n"
            laporan += f"Produk: {h['products']}\n"
            laporan += f"Energi Aktivasi: {h['Ea']}\n"
            laporan += f"Termodinamika: {h['thermo']}\n"
            laporan += "-"*40 + "\n"
        st.download_button("Unduh Laporan TXT",
                           data=laporan,
                           file_name="laporan_reaksi.txt",
                           mime="text/plain")
    else:
        st.info("Belum ada reaksi yang bisa disimpan.")
