# Clasificación de Cáncer de Mama: Naive Bayes vs k-NN

Comparación rigurosa de algoritmos de Machine Learning para diagnóstico 
automático mediante citología FNA del dataset Wisconsin Breast Cancer.

## 🏆 Resultado Principal
**k-NN con distancia Manhattan alcanza 98.04% de accuracy** 
(solo 4 errores en 204 casos), superando a Naive Bayes Multinomial (96.08%).

## 📊 Contenido
- Análisis exploratorio con t-SNE 3D
- Test χ² de independencia para selección de variables
- Comparación formal NB Multinomial vs k-NN
- Optimización de hiperparámetros (validación cruzada)
- Visualizaciones de fronteras de decisión con gradientes de confianza

## 🔬 Métodos
- Dataset: Wisconsin Breast Cancer (n=683, 9 variables)
- Split: 70% train / 30% test (n=204, evaluación independiente)
- Modelos: Naive Bayes Multinomial + k-NN (Euclidiana y Manhattan)

[Ver análisis completo en RPubs](tu-link-aqui)
