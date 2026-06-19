Las funciones Calibracion_General.cpp, General.cpp y Sort_GAGGs.cpp hay que incluirlas dentro de una carpeta src.
EL archivo GAGGs_Analysis.cpp es el que hay que ejecutar, y en él se selecciona la fuente utilizada y si se quiere realizar una calibración o un análisis de los datos ya calibrados.
El archivo General.cpp contiene las funciones utilizadas por el resto de funciones.
Calibracion_General.cpp realiza la calibración del espectro después de conocer la fuente medida. Se debe realizar antes de ejecutar Sort_GAGGs.cpp
Sort_GAGGs.cpp realiza el análisis del espectro, calculando la resolución, la eficiencia y realizando un análisis de add-back.
