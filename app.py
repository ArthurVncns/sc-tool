import streamlit as st

pg = st.navigation(
    {
        "": [
            st.Page("pages/home.py", title="Overview", default=True),
        ],
        "Analysis": [
            st.Page("pages/01_Upload.py", title="Upload"),
            st.Page("pages/02_QC.py", title="Quality Control"),
            st.Page("pages/03_Preprocessing.py", title="Preprocessing"),
            st.Page("pages/04_Reduction.py", title="Dimensionality Reduction"),
            st.Page("pages/05_Annotation.py", title="Annotation"),
            st.Page("pages/06_DE.py", title="Differential Expression"),
            st.Page("pages/07_Export.py", title="Export"),
        ],
    }
)
pg.run()
