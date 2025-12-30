# Clasificación de Cáncer de Mama: Naive Bayes Multinomial vs k-NN

[![R](https://img.shields.io/badge/R-4.5.2-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Status](https://img.shields.io/badge/Status-Complete-success.svg)]()

> **Comparación rigurosa de algoritmos de Machine Learning para diagnóstico automático de cáncer de mama mediante citología FNA del dataset Wisconsin Breast Cancer.**

---

## 🌐 **Ver Análisis Completo**

### 📊 [**ANÁLISIS INTERACTIVO COMPLETO →**](https://alej2andro.github.io/Clasificacion-Cancer-Mama-Naive-Bayes-Multinomial-vs-KNN/)

*Incluye visualizaciones 3D interactivas, fronteras de decisión con gradientes de confianza y análisis estadístico exhaustivo.*

---

## 🏆 **Resultado Principal**

**k-NN con distancia Manhattan alcanza 98.04% de accuracy** (solo 4 errores en 204 casos), superando significativamente a Naive Bayes Multinomial (96.08%).

| Modelo | Accuracy | Sensibilidad | Especificidad | Errores Totales |
|:-------|:--------:|:------------:|:-------------:|:---------------:|
| **k-NN Manhattan (k=17)** | **98.04%** | **97.18%** | **98.50%** | **4** |
| Naive Bayes Multinomial | 96.08% | 94.37% | 96.99% | 8 |

**Diferencia clave:** k-NN detecta **2 cánceres más** y genera **2 falsos positivos menos** que Naive Bayes.

---

## 📋 **Contenido del Análisis**

### 🔬 **Análisis Exploratorio**
- Visualización t-SNE 3D interactiva de separabilidad de clases
- Test χ² de Pearson para validación estadística de predictores
- Análisis de correlaciones y detección de outliers
- Matrices de confusión comparativas con mapas de calor

### 🤖 **Modelos Implementados**
1. **Naive Bayes Multinomial**
   - Configuración: 9 variables, `usekernel=FALSE`, `laplace=1`
   - Validación cruzada 10-fold
   - Análisis de probabilidades a posteriori
   - Fronteras de decisión en espacio 2D

2. **k-Nearest Neighbors**
   - Optimización de hiperparámetros (k=1 a k=35)
   - Comparación Euclidiana vs Manhattan
   - k óptimo: **17** (Manhattan)
   - Gradientes de confianza P(Maligno)

### 📊 **Visualizaciones Destacadas**
- Curvas ROC y Precision-Recall (AUC > 0.97)
- Radar chart multidimensional de métricas
- Fronteras de decisión con gradientes de confianza
- Mapa de aciertos/errores por caso individual
- Análisis de discrepancias entre algoritmos

---

## 🗂️ **Estructura del Dataset**

- **Fuente:** Wisconsin Breast Cancer Database (Dr. William H. Wolberg, 1992)
- **Observaciones:** 683 casos (444 benignos, 239 malignos)
- **Variables:** 9 características citológicas en escala ordinal 1-10
  - `Cl.thickness`, `Cell.size`, `Cell.shape`, `Marg.adhesion`
  - `Epith.c.size`, `Bare.nuclei`, `Bl.cromatin`
  - `Normal.nucleoli`, `Mitoses`
- **Split:** 70% train (n=479) / 30% test (n=204)

---

## 🛠️ **Metodología**

### **Preprocesamiento**
- Eliminación de valores faltantes (NA)
- Conversión a factores ordenados (escala 1-10)
- Estandarización para k-NN (center + scale)

### **Validación**
- Validación cruzada estratificada (10-fold)
- Test set completamente independiente (n=204)
- Seed fijada para reproducibilidad (`set.seed(123)`)

### **Métricas Evaluadas**
- Accuracy, Sensitivity, Specificity
- Precision (PPV), NPV, F1-Score
- Balanced Accuracy, Kappa
- Curvas ROC-AUC y PR-AUC
- Test de McNemar

---

## 💻 **Tecnologías y Paquetes**

### **Lenguaje**
- R 4.5.2

### **Paquetes Principales**
```r
# Datos y modelado
library(mlbench)      # Dataset Wisconsin Breast Cancer
library(caret)        # Train/test split, validación cruzada
library(naivebayes)   # Naive Bayes Multinomial
library(class)        # k-NN base
library(kknn)         # k-NN con distancia Manhattan

# Evaluación
library(pROC)         # Curvas ROC
library(PRROC)        # Curvas Precision-Recall

# Visualización
library(ggplot2)      # Gráficos base
library(plotly)       # t-SNE 3D interactivo
library(Rtsne)        # Reducción dimensional
library(corrplot)     # Matriz de correlaciones
library(fmsb)         # Radar chart
library(gridExtra)    # Composición de gráficos

# Manipulación
library(dplyr)        # Transformaciones
library(tidyr)        # Datos tidy
```

---

## 📂 **Estructura del Repositorio**
```
📦 naive-bayes-knn-ml-comparison/
├── 📄 Cancer_multinominal.Rmd          # Documento principal R Markdown
├── 📄 Cancer_multinominal.html         # Salida HTML compilada
├── 🖼️ Invasive_ductal_carcinoma.jpg   # Imagen histológica
├── 📄 README.md                         # Este archivo
└── 📄 LICENSE                           # Licencia MIT
```

---

## 🎯 **Hallazgos Clave**

### **1. Superioridad de k-NN Manhattan**
- **+1.96 pp** en Accuracy vs Naive Bayes
- **-2 falsos negativos** (detecta más cánceres)
- **-4 falsos positivos** (menos alarmas falsas)

### **2. Importancia de la Distancia Manhattan**
- Manhattan (98.04%) > Euclidiana (97.55%)
- Mayor robustez ante outliers citológicos
- Mejor especificidad (+0.76 pp)

### **3. Validación del Supuesto Naive**
- Correlación Cell.size ↔ Cell.shape: **r=0.907**
- A pesar de violar independencia, NB alcanza 96.08%
- k-NN libre de supuestos logra rendimiento superior

### **4. Variables Más Discriminantes**
1. `Cell.size` (χ²=540.02)
2. `Cell.shape` (χ²=523.21)
3. `Bare.nuclei` (χ²=489.32)

---

## 🔬 **Reproducibilidad**

### **Clonar y Ejecutar**
```bash
# Clonar repositorio
git clone https://github.com/alej2andro/naive-bayes-knn-ml-comparison.git
cd naive-bayes-knn-ml-comparison

# Abrir en RStudio y ejecutar
# File > Open File > Cancer_multinominal.Rmd
# Luego: Knit to HTML
```

### **Requisitos**
- R ≥ 4.0.0
- RStudio (recomendado)
- Paquetes listados arriba

---

## 📖 **Citar Este Trabajo**
```bibtex
@misc{figueroa2025breastcancer,
  author = {Figueroa Rojas, Alejandro},
  title = {Clasificación de Cáncer de Mama: Comparación Naive Bayes vs k-NN},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/alej2andro/naive-bayes-knn-ml-comparison}
}
```

---

## 📧 **Contacto**

**Autor:** Alejandro Figueroa Rojas  
**GitHub:** [@alej2andro](https://github.com/alej2andro)  
**Email:** [alejandro.figueroa.rojas@gmail.com](mailto:alejandro.figueroa.rojas@gmail.com)
**Proyecto:** [Ver Análisis Completo](https://alej2andro.github.io/Clasificacion-Cancer-Mama-Naive-Bayes-Multinomial-vs-KNN/)

---

## 🙏 **Referencias**

- **Dataset:** Wolberg, W. H. (1992). *Wisconsin Breast Cancer Database*. UCI Machine Learning Repository. [DOI: 10.24432/C5HP4Z](https://doi.org/10.24432/C5HP4Z)
- **Metodología:** Hastie, T., Tibshirani, R., & Friedman, J. (2009). *The Elements of Statistical Learning* (2nd ed.). Springer.

---

**⭐ Si este proyecto te resultó útil, considera darle una estrella en GitHub**
