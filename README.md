
# CSAR — X-Mapping Model for ²¹⁰Pb Dating

**Python Implementation — Release 1**  
Created by **J.M. Abril-Hernández (2025)**  
Departamento de Física Aplicada I, E.T.S.I. Agronómica, University of Seville (Spain)  
ORCID: [0000-0003-2540-5576](https://orcid.org/0000-0003-2540-5576)  
Email: [jmabril@us.es](mailto:jmabril@us.es)

---

## 📖 Description
Python implementation of the upgraded χ-mapping CSAR model for ²¹⁰Pb-based dating of recent sediments.  
Computes chronologies from empirical profiles using parametric domain search and confidence region analysis, with uncertainty estimates.

---

## 🔗 Related Publications
- Abril-Hernández, J.M. (2025). *Review on ²¹⁰Pb-based dating models for recent sediments*.  
  [https://doi.org/10.1016/j.jenvrad.2025.107749](https://doi.org/10.1016/j.jenvrad.2025.107749)  
- Abril, J.M. (2023). *²¹⁰Pb-based dating of recent sediments with the χ-mapping CSAR model*.  
  [https://doi.org/10.1016/j.jenvrad.2023.107247](https://doi.org/10.1016/j.jenvrad.2023.107247)

---

## ✅ Required Files
- `Core_C1.txt` — Empirical data with the ²¹⁰Pb profile.
- `/aleat_S1/` — Folder containing the library of random samples.
- `Random_generator_S1.py` — Script for generating the random sample library.
- `configuration.json` — Input file with all parameters.
- `CSAR_map.py` — Generates the 4D map for X_df.
- `CSAR_cronos.py` — Defines confidence region and outputs solutions.

---

## ▶️ How to Run
1. Prepare input data in `Core_C1.txt` (3-column format).
2. Generate random sample library:
   ```bash
   python Random_generator_S1.py
More details in README.PDF file within ZIP
