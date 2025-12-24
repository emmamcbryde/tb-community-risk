import streamlit as st
import os
import sys

# Ensure the project root is in sys.path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

# Now import your modules
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
