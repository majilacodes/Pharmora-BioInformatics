import streamlit as st
import json
import requests
from streamlit_lottie import st_lottie

# Set the page configuration
st.set_page_config(page_title="Landing Page", layout="wide")

# Center the title and subtitle using HTML with custom styles
st.markdown("""
    <h1 style='text-align: center; font-size: 8em; font-family: "Roboto", sans-serif; margin-top: -35px;'>Pharmora</h1>
    <h2 style='text-align: center; font-size: 1.9em; font-family: "Arial", sans-serif; margin-top: -45px;'>Accelerating Innovation in BioInformatics</h2>
""", unsafe_allow_html=True)

def load_lottieurl(url: str):
    r = requests.get(url)
    if r.status_code != 200:
        return None
    return r.json()

lottie_hello = load_lottieurl("https://lottie.host/5c7801eb-9d78-4186-aa86-5e988f0513af/O8O2lxLtNt.json")

st_lottie(
    lottie_hello,
    speed = 0.25,
    reverse = False,
    loop = True,
    quality = "high",
    height = "auto",
    width = "auto",
    key = None,
)

# Add an image
st.image("https://via.placeholder.com/800x400", caption="Aesthetic Image", use_column_width=True)

# Add a button
if st.button("Get Started"):
    st.write("Thank you for clicking the button!")

# Add a text input for user feedback
user_input = st.text_input("What do you think about this page?", "")

if user_input:
    st.write(f"Your feedback: {user_input}")

# Add a select box for choosing a theme
theme = st.selectbox("Choose a theme:", ["Light", "Dark", "Colorful"])

st.write(f"You selected the {theme} theme.")

# Add a slider for rating
rating = st.slider("Rate your experience (1-10)", 1, 10)

st.write(f"Your rating: {rating}")

# Add a sidebar for additional navigation or information
st.sidebar.header("Navigation")
st.sidebar.write("Use this sidebar to navigate through the app.")

# Create a section for cards
st.header("Explore Our Features")

# Create columns for cards
col1, col2, col3 = st.columns(3)

with col1:
    st.subheader("Feature 1")
    st.image("https://via.placeholder.com/150", caption="Feature 1 Image")
    st.write("Description of Feature 1.")
    if st.button("Learn More", key="feature1"):
        st.write("You clicked Learn More for Feature 1!")

with col2:
    st.subheader("Feature 2")
    st.image("https://via.placeholder.com/150", caption="Feature 2 Image")
    st.write("Description of Feature 2.")
    if st.button("Learn More", key="feature2"):
        st.write("You clicked Learn More for Feature 2!")

with col3:
    st.subheader("Feature 3")
    st.image("https://via.placeholder.com/150", caption="Feature 3 Image")
    st.write("Description of Feature 3.")
    if st.button("Learn More", key="feature3"):
        st.write("You clicked Learn More for Feature 3!")
