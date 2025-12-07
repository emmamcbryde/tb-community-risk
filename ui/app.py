import streamlit as st

from ui.dynamic_ui import render_dynamic_ui
from ui.static_ui import render_static_ui

st.title("🧮 TB Community Risk Models")

model_choice = st.radio(
    "Select model type:",
    ["Dynamic Model", "Static Model"]
)

if model_choice == "Dynamic Model":
    render_dynamic_ui()
else:
    render_static_ui()
