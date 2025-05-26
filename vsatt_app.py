import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

st.set_page_config(page_title="Vibration Analysis & Teaching Tool", layout="wide")

st.header("Vibration Analysis and Teaching Tool (VSATT)")

st.subheader("System Schematic")
st.image(
    "https://ctms.engin.umich.edu/CTMS/Content/Introduction/System/Modeling/figures/mass_spring_damper.png",
    caption="Fig. 1: Mass-Spring-Damper System"
)
st.image(
    "https://ctms.engin.umich.edu/CTMS/Content/Introduction/System/Modeling/figures/mass_spring_damper_FBD.png",
    caption="Fig. 2: Free Body Diagram of Mass-Spring-Damper System"
)

st.sidebar.header("Case Study 1 Parameters")
m1 = st.sidebar.slider('Mass 1 (kg)', 1.0, 100.0, 10.0, 1.0)
c1 = st.sidebar.slider('Damping 1 (Ns/m)', 0.0, 50.0, 4.0, 0.5)
k1 = st.sidebar.slider('Spring Constant 1 (N/m)', 10.0, 1000.0, 100.0, 10.0)
F1 = st.sidebar.slider('Input Force 1 (N)', 0.0, 500.0, 0.0, 10.0)
x01 = st.sidebar.slider('Initial Displacement 1 (m)', -1.0, 1.0, 0.01, 0.01)
v01 = st.sidebar.slider('Initial Velocity 1 (m/s)', -10.0, 10.0, 0.0, 0.1)

st.sidebar.header("Case Study 2 Parameters")
m2 = st.sidebar.slider('Mass 2 (kg)', 1.0, 100.0, 10.0, 1.0)
c2 = st.sidebar.slider('Damping 2 (Ns/m)', 0.0, 50.0, 10.0, 0.5)
k2 = st.sidebar.slider('Spring Constant 2 (N/m)', 10.0, 1000.0, 100.0, 10.0)
F2 = st.sidebar.slider('Input Force 2 (N)', 0.0, 500.0, 0.0, 10.0)
x02 = st.sidebar.slider('Initial Displacement 2 (m)', -1.0, 1.0, 0.01, 0.01)
v02 = st.sidebar.slider('Initial Velocity 2 (m/s)', -10.0, 10.0, 0.0, 0.1)

# Shared parameter
t_end = st.sidebar.slider('Simulation Time (s)', 1.0, 30.0, 10.0, 1.0)
t = np.linspace(0, t_end, 500)

# Function to compute response
def compute_response(m, c, k, F, x0, v0):
    wn = np.sqrt(k / m)
    zeta = c / (2 * np.sqrt(k * m))
    x0 += F / k if k != 0 else 0
    if zeta > 1:
        s1 = -wn * (zeta - np.sqrt(zeta ** 2 - 1))
        s2 = -wn * (zeta + np.sqrt(zeta ** 2 - 1))
        A = (v0 - s2 * x0) / (s1 - s2)
        B = x0 - A
        x = A * np.exp(s1 * t) + B * np.exp(s2 * t)
        label = f"Overdamped (ζ={zeta:.2f})"
        mode = "Overdamped"
    elif zeta == 1:
        A = x0
        B = v0 + wn * x0
        x = (A + B * t) * np.exp(-wn * t)
        label = f"Critically Damped (ζ={zeta:.2f})"
        mode = "Critically Damped"
    else:
        wd = wn * np.sqrt(1 - zeta ** 2)
        A = x0
        B = (v0 + zeta * wn * x0) / wd
        x = np.exp(-zeta * wn * t) * (A * np.cos(wd * t) + B * np.sin(wd * t))
        label = f"Underdamped (ζ={zeta:.2f})"
        mode = "Underdamped"
    return x, label, wn, zeta, mode

x1, label1, wn1, zeta1, mode1 = compute_response(m1, c1, k1, F1, x01, v01)
x2, label2, wn2, zeta2, mode2 = compute_response(m2, c2, k2, F2, x02, v02)

# Plot both case studies
st.subheader("Displacement Response Comparison")
fig, ax = plt.subplots()
ax.plot(t, x1, label=f"Case 1 - {label1}")
ax.plot(t, x2, label=f"Case 2 - {label2}", linestyle='--')
ax.set_xlabel("Time (s)")
ax.set_ylabel("Displacement (m)")
ax.set_title("Comparison of Case Study 1 and 2")
ax.grid(True)
ax.legend()
st.pyplot(fig)

# Show natural frequency and damping ratio
st.subheader("Summary of Case Studies")
st.markdown(f"**Case Study 1:** ωₙ = {wn1:.3f} rad/s, ζ = {zeta1:.3f}, Mode = {mode1}")
st.markdown(f"**Case Study 2:** ωₙ = {wn2:.3f} rad/s, ζ = {zeta2:.3f}, Mode = {mode2}")

st.subheader("Reference")
st.write("Adapted from standard vibration system theories. Developed for teaching use.")

st.subheader("Assessment Tools")
st.markdown('[Pre-Test](https://forms.gle/wPvcgnZAC57MkCxN8)', unsafe_allow_html=True)
st.markdown('[Post-Test](https://forms.gle/FdiKqpMLzw9ENscA9)', unsafe_allow_html=True)
