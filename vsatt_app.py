import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

st.set_page_config(page_title="VSATT - Vibration Teaching Tool", layout="wide")

# Fungsi simulasi getaran
def vibration_response(mass, damping, spring, force, x0, v0, t_end):
    k = spring
    c = damping
    m = mass
    wn = np.sqrt(k / m)
    zeta = c / (2 * np.sqrt(k * m))
    t = np.linspace(0, t_end, 500)

    if zeta > 1:
        s1 = -wn * (zeta - np.sqrt(zeta ** 2 - 1))
        s2 = -wn * (zeta + np.sqrt(zeta ** 2 - 1))
        A = (v0 - s2 * x0) / (s1 - s2)
        B = x0 - A
        x = A * np.exp(s1 * t) + B * np.exp(s2 * t)
        mode = "Overdamped"
    elif zeta == 1:
        A = x0
        B = v0 + wn * x0
        x = (A + B * t) * np.exp(-wn * t)
        mode = "Critically Damped"
    else:
        wd = wn * np.sqrt(1 - zeta ** 2)
        A = x0
        B = (v0 + zeta * wn * x0) / wd
        x = np.exp(-zeta * wn * t) * (A * np.cos(wd * t) + B * np.sin(wd * t))
        mode = "Underdamped"

    return t, x, wn, zeta, mode

st.title("VSATT - Vibration Simulation & Analysis Teaching Tool")

st.subheader("Input Parameters")

col1, col2 = st.columns(2)

inputs = []

for i, col in enumerate([col1, col2]):
    with col:
        st.markdown(f"### Case Study {i + 1}")
        mass = st.number_input(f"Mass {i + 1} (kg)", value=10.0, key=f"mass{i}")
        damping = st.number_input(f"Damping Constant {i + 1} (Ns/m)", value=10.0 if i == 0 else 4.0, key=f"damp{i}")
        spring = st.number_input(f"Spring Constant {i + 1} (N/m)", value=1.0, key=f"spring{i}")
        force = st.number_input(f"Force {i + 1} (N)", value=1.0, key=f"force{i}")
        x0 = st.number_input(f"Initial Position {i + 1}", value=1.0, key=f"x0{i}")
        v0 = st.number_input(f"Initial Velocity {i + 1}", value=1.0, key=f"v0{i}")
        t_end = st.number_input(f"Simulation Time {i + 1} (s)", value=100.0, key=f"t{i}")

        inputs.append((mass, damping, spring, force, x0, v0, t_end))

if st.button("Simulate Both Cases"):
    fig, ax = plt.subplots(figsize=(10, 5))
    
    for i, params in enumerate(inputs):
        t, x, wn, zeta, mode = vibration_response(*params)
        ax.plot(t, x, label=f"Case {i+1} - {mode}")
        st.markdown(f"### Output for Case Study {i + 1}")
        st.markdown(f"**Natural Frequency (wn):** {wn:.3f} rad/s")
        st.markdown(f"**Damping Ratio (zeta):** {zeta:.3f}")
        st.markdown(f"**System Type:** {mode}")

    ax.set_title("Comparison of Vibration Response")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Displacement (m)")
    ax.grid(True)
    ax.legend()
    st.pyplot(fig)
