#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 10:56:57 2023

@author: usuario
"""



import os
import shutil

def copy_and_rename_files(source_directory, destination_directory):
    # Verifica si el directorio de destino existe, si no, lo crea
    if not os.path.exists(destination_directory):
        os.makedirs(destination_directory)

    # Recorre por todos los directorios, subdirectorios y archivos en el directorio fuente
    for foldername, subfolders, filenames in os.walk(source_directory):
        for filename in filenames:
            if filename.startswith("short_summary"):
                # Construye la ruta completa del archivo original y de su destino
                source_path = os.path.join(foldername, filename)
                dest_filename = f"short_summary_{os.path.basename(foldername)}"
                dest_path = os.path.join(destination_directory, dest_filename)
                
                # Copia y renombra el archivo al directorio de destino
                shutil.copy2(source_path, dest_path)

src_dir = "/home/usuario/Documentos/Carlos_PhD/Nocardiopsis/BUSCO_streptosporangiales"  # Reemplaza esto con la ruta a tu directorio principal
dest_dir = "/home/usuario/Documentos/Carlos_PhD/Nocardiopsis/BUSCO_streptosporangiales/short_summaries"   # Reemplaza esto con la ruta a donde quieres guardar los archivos copiados

copy_and_rename_files(src_dir, dest_dir)



import os
import re

def extract_values_from_line(line):
    # Expresión regular para extraer los valores
    pattern = r"C:(?P<C>[\d\.]+)%\[S:(?P<S>[\d\.]+)%,D:(?P<D>[\d\.]+)%\],F:(?P<F>[\d\.]+)%,M:(?P<M>[\d\.]+)%"
    match = re.search(pattern, line)

    if match:
        return match.groupdict()  # Devuelve un diccionario con los valores extraídos
    return None

def extract_values_from_files(directory):
    results = []

    # Itera sobre cada archivo en el directorio
    for filename in os.listdir(directory):
        if filename.startswith("short_summary"):
            filepath = os.path.join(directory, filename)

            with open(filepath, 'r') as file:
                for current_line_number, line in enumerate(file, 1):
                    if current_line_number == 9:
                        values = extract_values_from_line(line)
                        if values:
                            # Agrega el nombre del archivo al diccionario
                            values["filename"] = filename
                            results.append(values)
                        break

    return results


directory_path = "/home/usuario/Documentos/Carlos_PhD/Nocardiopsis/BUSCO_streptosporangiales/short_summaries"  # Reemplaza esto con la ruta donde guardaste los archivos copiados
values_list = extract_values_from_files(directory_path)
    
# Crear listas vacías para cada clave
C_values = []
S_values = []
D_values = []
F_values = []
M_values = []

# Iterar sobre cada diccionario en values_list y extraer los valores
for item in values_list:
    C_values.append(float(item["C"]))
    S_values.append(float(item["S"]))
    D_values.append(float(item["D"]))
    F_values.append(float(item["F"]))
    M_values.append(float(item["M"]))

# Imprimir los valores (puedes omitir esto si no es necesario)
print("Valores de C:", C_values)
print("Valores de S:", S_values)
print("Valores de D:", D_values)
print("Valores de F:", F_values)
print("Valores de M:", M_values)


import matplotlib.pyplot as plt
import seaborn as sns

# Asumiendo que ya has extraído los valores de C en la lista C_values...

# Configurar el estilo de Seaborn
sns.set_style("whitegrid")

# Crear una figura con dos subplots: uno para el histograma y otro para el boxplot
fig, ax = plt.subplots(nrows=2, figsize=(10, 8))

# Histograma
ax[0].hist(C_values, bins=20, edgecolor='black', color='blue', alpha=0.7)
ax[0].set_title('Histograma de C')
ax[0].set_xlabel('Valor de C')
ax[0].set_ylabel('Frecuencia')

# Boxplot
sns.boxplot(x=C_values, ax=ax[1], color='lightblue')
ax[1].set_title('Boxplot de C')
ax[1].set_xlabel('Valor de C')

# Ajustar el layout y mostrar la figura
plt.tight_layout()


#outliers_files
outliers_files = [item["filename"] for item in values_list if float(item["C"]) < 95]
outliers_files2 = [item["filename"].replace("short_summary_", "") for item in values_list if float(item["C"]) < 95]

import os

# Define la ruta de la carpeta donde se encuentran los archivos a eliminar
folder_path = "/home/usuario/Documentos/Carlos_PhD/Nocardiopsis/raw_data/Nocardiopsis_fasta_BUSCO"

for filename in outliers_files2:
    # Añade el prefijo "short_summary_" al nombre del archivo
    full_filename = filename+".fna"
    
    # Construye la ruta completa al archivo
    file_path = os.path.join(folder_path, full_filename)

    # Verifica si el archivo existe antes de intentar eliminarlo
    if os.path.exists(file_path):
        os.remove(file_path)
        print(f"Archivo {full_filename} eliminado con éxito.")
    else:
        print(f"El archivo {full_filename} no se encontró.")
