## Intro
Disponemos de 2 modelos murinos de melanoma metastático cerebral en los que tenemos por un lado imágenes de cortes histológicos teñidas por DAPI (nucleos) y RNA-FISH de 100 transcritos, que permite la localización de marcadores de distintos tipos celulares y fenotipos a nivel de célula individual en el tejido.

Queremos conseguir una segmentación adecuada o hacer agrupaciones con los distintos tipos celulares basándonos principalmente en los 100 marcadores celulares que tenemos. 

### Problema:  

Hay dos tipos de segmentaciones que se pueden realizar a través del software de RESOLVE. Sin embargo ninguna está cerca de una segmentación aceptable.  

- La segmentación por DAPI genera, para todas las células, un área alrededor del núcleo basada en micras. Resultando en que todas las células tienen siempre el mismo tamaño de área alrededor y esto no es correcto.  

- La segmentación por señales de RNA genera células sin considerar el núcleo de éstas, que en casi todos los casos se queda fuera de la célula porque las sondas que tenemos marcadas están en el citoplasma ya que pertenecen a RNA de transcritos. 

Además, para ambos tipos de segmentación las zonas con mucho solapamiento no hay forma de predecirlas. 

#### Otros problemas que surgieron de ciertas muestras:
- El área de contaminación en la tile 13 de la muestra B1-2, con expresión de Cd8. 
- El infiltrado generalizado de Cd4 en la muestra C2-1. Que parece ser debido a expresión de Cd4 por células neuronales.

## Desarrollo Agglomerative Clustering 
Vamos a usar el algoritmo de sklearn [`aggomerativeClustering`](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.AgglomerativeClustering.html). Introduciremos las coordenadas de los transcritos, para ello tenemos que translocar los datos para conseguir una matriz de expresión binaria por coordenadas, donde cada fila es un punto con coordenadas $x$ e $y$, y cada columna indica si un gen se expresa (1) o no (0) en ese punto. Vamos a juntar las coordenadas z en un mismo plano, y por lo tanto juntaremos las que coincidan, para reducir información.

### Codigo para reformateo de datos

```py
from collections import defaultdict

input_file = "C:\\Users\\user\\Escritorio\\practicas\\32873-slide3_B1-2_results.txt"
output_file = "C:\\Users\\user\\Escritorio\\practicas\\procesadosB1-2.csv"

# Dictionary with keys = (x, y), values = set de genes expresados
data = defaultdict(set)
all_genes = set()

# Leer archivo y llenar estructuras
with open(input_file, 'r') as f_in:
    for line in f_in:                     
        columns = line.strip().split()
        if len(columns) < 4:
            continue  # Jump imcomplete lines (just in case)
        x = columns[0]
        y = columns[1]
        gen = columns[3]
        all_genes.add(gen)
        data[(x, y)].add(gen)

# Order genes for consistency in columns
gene_list = sorted(all_genes)


print(f"Total coordenadas únicas: {len(data)}")
print(f"Total genes únicos: {len(gene_list)}")
print(gene_list)
```

#### Resultados
```
Total coordenadas únicas: 1698725
Total genes únicos: 100
['Acta2', 'Adgre1', 'Aldh1l1', 'Aldoc', 'Aqp4', 'Arg1', 'Axl', 'Ccl2', 'Ccl3', 'Ccl5', 'Ccr2', 'Cd19', 'Cd24a', 'Cd274', 'Cd3e', 'Cd4', 'Cd44', 'Cd47', 'Cd52', 'Cd69', 'Cd74', 'Cd8a', 'Csf1', 'Csf1r', 'Csf3r', 'Cspg4', 'Ctla4', 'Cxcl10', 'Cxcl12', 'Cxcl16', 'Cxcl3', 'Cxcr3', 'Cxcr4', 'Dct', 'Flt3', 'Fn1', 'Foxp3', 'Gap43', 'Gfap', 'Grem1', 'Gzmb', 'H2-K1', 'Havcr2', 'Hmga2', 'Icam1', 'Ifng', 'Ifngr1', 'Il10', 'Il1b', 'Itga4', 'Itgae', 'Itgam', 'Itgax', 'Itgb1', 'Klrg1', 'Lag3', 'Ly6g', 'Map2', 'Mitf', 'Mki67', 'Mmp9', 'Mog', 'Mrc1', 'Myh11', 'Ncr1', 'Nefm', 'Ngfr', 'Nkx1-1', 'Nr2f2', 'Nsg2', 'Olig1', 'P2ry12', 'Pdcd1', 'Pecam1', 'Plxnc1', 'Pmel', 'Prf1', 'Ptprc', 'Rbfox3', 'Sema7a', 'Serpine1', 'Sirpa', 'Slc7a5', 'Sox10', 'Sox4', 'Sox9', 'Sparcl1', 'Spp1', 'Stat1', 'Stat3', 'Tbx2', 'Tcf4', 'Tgfb2', 'Tgfbr1', 'Timp1', 'Tlr4', 'Tmem119', 'Trem2', 'Vegfa', 'Zeb1']
```

#### Output file
```py
# Output file
with open(output_file, 'w') as f_out:

#csv
    # header
    f_out.write("x,y," + ",".join(gene_list) + "\n") # csv
    for (x, y), expressed_genes in data.items():
        row = [x, y]
        row += ['1' if g in expressed_genes else '0' for g in gene_list]
        f_out.write(",".join(row) + "\n")
```

#### Visualizacion

```py
import pandas as pd
import matplotlib.pyplot as plt

output_file = "C:\\Users\\user\\Escritorio\\practicas\\procesadosB1-2.csv"

# Cargar archivo generado
file1 = pd.read_csv(output_file)

plt.figure(figsize=(9, 15))
plt.scatter(file1['x'], file1['y'], s=1, alpha=0.2)
plt.xlabel("Coordenada X")
plt.ylabel("Coordenada Y")
plt.title("Distribución de todos los genes")
plt.grid(True)
plt.tight_layout()
plt.show()
```

<img width="889" height="1490" alt="image" src="https://github.com/user-attachments/assets/fbd8f43b-4360-4a26-9518-077d6a00cca0" />

*La imagen está formada por 6375 de ancho x 10625 de alto. Cada cuadricula son 2125 x 2125 pixeles.*

### TILES
Al intentar correr el algoritmo de clustering con esta imagen, esto era demasiado pesado y la memoria no lo permitía, con lo que nos vimos en la situacion de reducirla a las tiles que componen la imagen.

#### Estrategia para obtener las tiles

1. Leer el archivo procesado (procesadosB1-2.csv).
2. Obtener los valores mínimos y máximos de x e y.
3. Dividir el rango de x en 3 partes (tiles en x) y el de y en 5 partes (tiles en y).
4. Asignar cada coordenada a su tile correspondiente.
5. Guardar un archivo por tile.

```py
import pandas as pd
import csv
import os
import matplotlib.pyplot as plt

input_file = "C:\\Users\\user\\Escritorio\\practicas\\procesadosB1-2.csv"
output_dir = "C:\\Users\\user\\Escritorio\\practicas\\tiles_B1-2"

os.makedirs(output_dir, exist_ok=True)

# === LECTURA ===
with open(input_file, newline='') as f:
    reader = list(csv.reader(f))
    header = reader[0]
    rows = reader[1:]

# Convertir coordenadas a enteros
for row in rows:
    row[0] = int(row[0])  # x
    row[1] = int(row[1])  # y

# Obtener rangos
x_vals = [row[0] for row in rows]
y_vals = [row[1] for row in rows]
x_min, x_max = min(x_vals), max(x_vals)
y_min, y_max = min(y_vals), max(y_vals)

# División en tiles
num_tiles_x = 3  # columnas
num_tiles_y = 5  # filas
tile_width = (x_max - x_min) / num_tiles_x
tile_height = (y_max - y_min) / num_tiles_y

# Guardar datos por tile
tiles = {}
resumen = {}

for row in rows:
    x, y = row[0], row[1]
    col = int((x - x_min) / tile_width) + 1
    row_inv = int((y - y_min) / tile_height)
    row_num = num_tiles_y - row_inv  # invertimos para contar desde arriba
    if col > num_tiles_x: col = num_tiles_x
    if row_num < 1: row_num = 1
    if row_num > num_tiles_y: row_num = num_tiles_y
    tile_id = (row_num, col)
    tiles.setdefault(tile_id, []).append(row)


# === GUARDAR TILES Y MOSTRAR RESUMEN ===
for (fila, col), tile_rows in sorted(tiles.items()):
    tile_filename = f"tile_{fila}_{col}.csv"
    tile_path = os.path.join(output_dir, tile_filename)

    with open(tile_path, 'w', newline='') as f_out:
        writer = csv.writer(f_out)
        writer.writerow(header)
        for row in tile_rows:
            writer.writerow([str(row[0]), str(row[1])] + row[2:])

    # Mostrar resumen
    genes_expresados = set()
    for row in tile_rows:
        genes_expresados.update([
            gene for gene, val in zip(header[2:], row[2:]) if val == "1"
        ])
    print(f"Tile {fila}_{col}: {len(tile_rows)} coordenadas, {len(genes_expresados)} genes expresados")

print(f"\nTiles guardados en: {output_dir}")
```

#### Resultados
```
Tile 1_1: 104915 coordenadas, 97 genes expresados
Tile 1_2: 125673 coordenadas, 94 genes expresados
Tile 1_3: 166618 coordenadas, 95 genes expresados
Tile 2_1: 81291 coordenadas, 88 genes expresados
Tile 2_2: 122005 coordenadas, 98 genes expresados
Tile 2_3: 187555 coordenadas, 95 genes expresados
Tile 3_1: 85039 coordenadas, 88 genes expresados
Tile 3_2: 103449 coordenadas, 96 genes expresados
Tile 3_3: 132953 coordenadas, 95 genes expresados
Tile 4_1: 69330 coordenadas, 89 genes expresados
Tile 4_2: 100256 coordenadas, 96 genes expresados
Tile 4_3: 113823 coordenadas, 96 genes expresados
Tile 5_1: 68973 coordenadas, 88 genes expresados
Tile 5_2: 121559 coordenadas, 96 genes expresados
Tile 5_3: 115286 coordenadas, 96 genes expresados

Tiles guardados en: C:\Users\user\Escritorio\practicas\tiles_B1-2
```

#### Visualizacion de la tile
Elegimos la tile con más variedad celular.

```py
import pandas as pd
import matplotlib.pyplot as plt

tile_id = (2, 2)
subtile_file = f"C:\\Users\\user\\Escritorio\\practicas\\tile_{tile_id[0]}_{tile_id[1]}.csv"

# Cargar archivo generado
file1 = pd.read_csv(subtile_file)
plt.figure(figsize=(5, 5))
plt.scatter(file1['x'], file1['y'], s=1, alpha=0.2)
plt.xlabel("Coordenada X")
plt.ylabel("Coordenada Y")
plt.title("Distribución de todos los genes")
plt.grid(True)
plt.tight_layout()
plt.show()
```

<img width="500" height="500" alt="image" src="https://github.com/user-attachments/assets/2ca441f9-3a9b-40e5-bf5d-b5d595052e53" />

### Subtiles
La imagen sigue siendo demasiado pesada para trabajar con ella. Volvemos a particionar, esta vez cada tile anterior la partimos en 4 de la siguiente manera:

<img width="434" height="619" alt="image" src="https://github.com/user-attachments/assets/9e266ae6-c435-4759-9eec-5f7891b285ed" />

```py
# === INPUT ===
tile_id = (2, 2)
tile_file = f"C:\\Users\\user\\Escritorio\\practicas\\tiles_B1-2\\tile_{tile_id[0]}_{tile_id[1]}.csv"
output_subdir = f"C:\\Users\user\\Escritorio\\practicas\\tiles_B1-2\\tile_{tile_id[0]}_{tile_id[1]}_subtiles"
os.makedirs(output_subdir, exist_ok=True)

# === LECTURA DE TILE ===
with open(tile_file, newline='') as f:
    reader = list(csv.reader(f))
    header = reader[0]
    rows = reader[1:]

# Coordenadas
for row in rows:
    row[0] = int(row[0])  # x
    row[1] = int(row[1])  # y

x_vals = [row[0] for row in rows]
y_vals = [row[1] for row in rows]
x_min, x_max = min(x_vals), max(x_vals)
y_min, y_max = min(y_vals), max(y_vals)

# Dimensiones
mid_x = (x_max + x_min) / 2
mid_y = (y_max + y_min) / 2

# Clasificación en subtiles
subtiles = {
    'A': [],  # top-left
    'B': [],  # top-right
    'C': [],  # bottom-left
    'D': []   # bottom-right
}

for row in rows:
    x, y = row[0], row[1]
    if y >= mid_y and x < mid_x:
        subtiles['A'].append(row)
    elif y >= mid_y and x >= mid_x:
        subtiles['B'].append(row)
    elif y < mid_y and x < mid_x:
        subtiles['C'].append(row)
    elif y < mid_y and x >= mid_x:
        subtiles['D'].append(row)

# Guardar subtiles
for label, subtile_rows in subtiles.items():
    out_path = os.path.join(output_subdir, f"tile_{tile_id[0]}_{tile_id[1]}_{label}.csv")
    with open(out_path, 'w', newline='') as f_out:
        writer = csv.writer(f_out)
        writer.writerow(header)
        for row in subtile_rows:
            writer.writerow([str(row[0]), str(row[1])] + row[2:])
    print(f"Subtile {label}: {len(subtile_rows)} coordenadas guardadas en {out_path}")
```

#### Resultado

```py
import pandas as pd
import matplotlib.pyplot as plt

tile_id = (2, 2, "D")
subtile_file = f"C:\\Users\\user\\Escritorio\\practicas\\tile_{tile_id[0]}_{tile_id[1]}_subtiles\\tile_{tile_id[0]}_{tile_id[1]}_{tile_id[2]}.csv"

# Cargar archivo generado
file1 = pd.read_csv(subtile_file)
plt.figure(figsize=(5, 5))
plt.scatter(file1['x'], file1['y'], s=1, alpha=0.2)
plt.xlabel("Coordenada X")
plt.ylabel("Coordenada Y")
plt.title("Distribución de todos los genes")
plt.grid(True)
plt.tight_layout()
plt.show()
```

<img width="490" height="490" alt="image" src="https://github.com/user-attachments/assets/6ffa6ba5-d74c-4b78-846f-c605b4d3a7a6" />

## Clustering

### Gene selection + Sample selection
```py
# Genes a usar 
genes_interes = ['Gfap','Nefm','Map2','Aldh1l1','Aldoc','Aqp4','Foxp3','Cd8a', 'Sparcl1', 
                 'Flt3','Itgae','Itgax','Ptprc','Dct','Mitf','Pmel','Slc7a5','P2ry12','Cd3e','Cd4', 
                 'Tmem119','Cd24a','Itgam','Ly6g','Adgre1','Nsg2','Rbfox3','Trem2','Klrg1', 
                 'Nkx1-1','Mog','Olig1','Acta2','Fn1','Myh11','Pecam1']  

# Tile a usar 
tile_id = (2, 2, "D") 
subtile_file = f"C:\\Users\\user\\Escritorio\\practicas\\tile_{tile_id[0]}_{tile_id[1]}_subtiles\\tile_{tile_id[0]}_{tile_id[1]}_{tile_id[2]}.csv" 

# Cargar archivo generado 
df_expr = pd.read_csv(subtile_file) 

# Asegurarse de que los genes existan en el archivo 
genes_existentes = [g for g in genes_interes if g in df_expr.columns] 
if not genes_existentes: 
    raise ValueError("Ninguno de los genes indicados está presente en el archivo.") 
 
# Filtrar filas donde al menos uno de los genes está expresado 
mask = df_expr[genes_existentes].astype(int).sum(axis=1) > 0 
df_filtrado = df_expr[mask].copy() 
```

### Colores para genes clave

```py
genes_coloreados = {
    # Astrocytes (azules)
    'Gfap': '#1f77b4',         # azul oscuro
    'Aldh1l1': '#4F9ED2',      # azul medio
    'Aldoc': '#77B7E5',        # azul claro
    'Aqp4': '#ADD8E6',         # azul muy claro
    'Sparcl1': '#5DADE2',      # azul vibrante

    # Neurons / Melanoma (marrón claro)
    'Map2': '#8e44ad',         # morado oscuro (mantener)
    'Nefm': '#9b59b6',         # morado medio (mantener)
    'Nsg2': '#BB8FCE',         # morado claro (mantener)
    'Rbfox3': '#C39BD3',       # morado claro (mantener)
    'Trem2': '#D7BDE2',        # morado muy claro (mantener)

    # Melanoma (marrones claros)
    'Dct': 'saddlebrown',   
    'Mitf': 'chocolate', 
    'Pmel': 'peru',  
    'Slc7a5': 'sandybrown',

    # T cells (rojos)
    'Cd4': '#B71C1C',
    'Foxp3': '#C62828',
    'Cd8a': '#D32F2F',
    'Cd3e': '#E53935',

    # B cells (rojos)
    'Cd19': '#F44336',

    # NK cells (rojos)
    'Klrg1': '#E53935',
    'Nkx1-1': '#F44336',

    # Leukocytes (rojos)
    'Ptprc': '#EF5350',

    # Dendritic cells (verde oscuro)
    'Flt3': '#145A32',
    'Itgae': '#1E8449',
    'Itgax': '#27AE60',

    # Microglia (verde oliva)
    'P2ry12': '#7D6608',
    'Tmem119': '#B7950B',

    # Mieloides (amarillos)
    'Cd24a': '#FBC02D',        # amarillo fuerte
    'Itgam': '#FDD835',        # amarillo medio
    'Ly6g': '#FFF176',         # amarillo suave
    'Adgre1': '#FFF9C4',       # amarillo muy pálido

    # Oligodendrocytes (cian)
    'Mog': '#17A589',
    'Olig1': '#48C9B0',

    # Vasculature (negro y grises)
    'Acta2': '#212121',
    'Fn1': '#616161',
    'Myh11': '#9E9E9E',
    'Pecam1': '#BDBDBD'
}
```

### Funcion para clustering

```py
def cluster_and_plot_cells_with_genes(df, batch, dataset, genes_coloreados, 
                                      distance_threshold=60, alpha_val=0.05):
    """
    Clustering espacial + overlay de genes coloreados + alpha shapes por cluster
    """

    # Filtrar datos
    subset = df[(df['batch_nr'] == batch) & (df['Dataset'] == dataset)].copy()
    if subset.empty:
        print(f"No data found for batch '{batch}' and dataset '{dataset}'.")
        return None

    # Clustering jerárquico
    coords = subset[['Coordenada X', 'Coordenada Y']]
    
    # Conectividad basada en vecinos más cercanos
    clustering = AgglomerativeClustering(
        n_clusters=None,
        distance_threshold=distance_threshold,
        linkage='average'
    )
    labels = clustering.fit_predict(coords)
    subset['cluster'] = labels

    print(f"Number of cells found: {len(set(labels))}")

    # === PLOTEO ===
    fig, ax = plt.subplots(figsize=(8, 8))

    # 1) Todos los puntos en gris
    ax.scatter(
        subset['Coordenada X'],
        subset['Coordenada Y'],
        s=1
    )

    # 2) Genes coloreados
    for gene, color in genes_coloreados.items():
        if gene not in df.columns:
            print(f"[AVISO] Gen {gene} no está en el DataFrame.")
            continue
        subset_gene = subset[subset[gene] == 1]
        ax.scatter(
            subset_gene['Coordenada X'],
            subset_gene['Coordenada Y'],
            s=1,
            color=color,
            label=gene
        )

    # 3) Dibujar contornos de los clusters
    for clust_id in subset['cluster'].unique():
        clust_data = subset[subset['cluster'] == clust_id]
        points = clust_data[['Coordenada X', 'Coordenada Y']].values

        if len(points) < 4:
            continue

        try:
            shape = alphashape.alphashape(points, alpha_val)
            shapes = [shape] if isinstance(shape, Polygon) else shape.geoms

            for geom in shapes:
                x, y = geom.exterior.xy
                ax.plot(x, y, color='black', linewidth=1)
        except Exception as e:
            print(f"Alpha shape failed for cluster {clust_id}: {e}")

    # 4) Estetica de visualización
    # Guardar como PNG sin bordes ni leyenda, con fondo transparente
    ax.set_title('')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks([])
    ax.set_yticks([])

    # Eliminar bordes (spines)
    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.savefig(
        f"cluster_plot_tile_{tile_id[0]}_{tile_id[1]}_{tile_id[2]}.png",
        dpi=300,
        transparent=True,
        bbox_inches='tight',
        pad_inches=0
    )

    plt.show()

    return subset


clustered_df = cluster_and_plot_cells_with_genes( 
    df_filtrado, 
    batch='B1', 
    dataset='slide1', 
    genes_coloreados=genes_coloreados, 
    distance_threshold=75 
) 
```

<img width="640" height="636" alt="image" src="https://github.com/user-attachments/assets/55ff716a-81e9-4ec5-bcfe-d2cf2c05ca2e" />

> Los resultados no son óptimos, se mezclan distintos genes en el mismo cluster y los tamaños no concuerdan con lo que esperaríamos.
> 
> Lo que ocurre es que hasta ahora estamos agrupando solo teniendo en cuenta la distribución espacial de las coordenadas, no los genes persé. Los cuales representan los tipos celulares. Necesitamos introducir la agrupación génica.
> 
> Probamos con una matriz de conectividad sin exito. Vamos a probar a normalizar y añadir clusterización génica a la espacial que hemos realizado hasta ahora.

## Añadiendo la Agrupación Espacial + Génica

Hasta ahora hemos hecho un clustering espacial, basado en coordenadas unicamente.

```py
coords = subset[['Coordenada X', 'Coordenada Y']]
```

Para ésto no necesitabamos normalizar, porque las coordenadas espaciales (x,y) representan posiciones reales en el tejido y su escala es significativa. Si las normalizamos, distorsionamos las distancias físicas entre células. 

Pero para añadir la clusterización génica hay que normalizar los datos.

```py
features = df[['x', 'y'] + genes] 
```

Ya que aquí tenemos dos tipos de variables: 

- Coordenadas en rangos de [0–2000] 
- Genes (binarios)

Con lo que las coordenadas dominan las distancias numéricas en el clustering. El clustering será casi 100% guiado por la posición, porque las diferencias en x, y son mucho mayores. 

Por eso en este caso sí debemos normalizar todo, tanto $x$ e $y$, como los genes. 

Entonces la "distancia" entre dos células será una combinación de: 
- Qué tan cerca están físicamente 
- Cuántos genes comparten 

Y el objetivo es que ninguna variable domine la distancia total. 

Tras varios intentos con distintos parámetros no se obtenía nada óptimo. 

Probamos a dar color cada grupo, con el color del gen dominante en éste y añadimos una ponderación a los genes escalados para que tengan mayor peso de cara al clustering. Jugando con los distintos parámetros finalmente obtuvimos algo con lo que creemos que podemos avanzar a los análisis de vecinos y otros:

### Reajuste de colores para ver linfocitos y otras células más claramente
```py
genes_coloreados = {
    # Astrocytes (azules)
    'Gfap': '#1f77b4',         # azul oscuro
    'Aldh1l1': '#4F9ED2',      # azul medio
    'Aldoc': '#77B7E5',        # azul claro
    'Aqp4': '#ADD8E6',         # azul muy claro
    'Sparcl1': '#5DADE2',      # azul vibrante

    # Neurons / Melanoma (marrón claro)
    'Map2': '#8e44ad',         # morado oscuro (mantener)
    'Nefm': '#9b59b6',         # morado medio (mantener)
    'Nsg2': '#BB8FCE',         # morado claro (mantener)
    'Rbfox3': '#C39BD3',       # morado claro (mantener)
    'Trem2': '#D7BDE2',        # morado muy claro (mantener)

    # Melanoma (marrones claros)
    'Dct': 'saddlebrown',   
    'Mitf': 'chocolate', 
    'Pmel': 'peru',  
    'Slc7a5': 'sandybrown',

    # T cells (rojos)
    'Cd4': '#B71C1C',
    'Foxp3': '#C62828',
    'Cd8a': '#D32F2F',
    'Cd3e': '#E53935',

    # B cells (rojos)
    'Cd19': '#F44336',

    # NK cells (rojos)
    'Klrg1': '#E53935',
    'Nkx1-1': '#F44336',

    # Leukocytes (rojos)
    'Ptprc': '#EF5350',

    # Dendritic cells (verde oscuro)
    'Flt3': '#145A32',
    'Itgae': '#1E8449',
    'Itgax': '#27AE60',

    # Microglia (verde oliva)
    'P2ry12': '#7D6608',
    'Tmem119': '#B7950B',

    # Mieloides (amarillos)
    'Cd24a': '#FBC02D',        # amarillo fuerte
    'Itgam': '#FDD835',        # amarillo medio
    'Ly6g': '#FFF176',         # amarillo suave
    'Adgre1': '#FFF9C4',       # amarillo muy pálido

    # Oligodendrocytes (cian)
    'Mog': '#17A589',
    'Olig1': '#48C9B0',

    # Vasculature (negro y grises)
    'Acta2': '#212121',
    'Fn1': '#616161',
    'Myh11': '#9E9E9E',
    'Pecam1': '#BDBDBD'
}
```

### Librerias necesarias
```py
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering
from sklearn.neighbors import radius_neighbors_graph
import alphashape
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import numpy as np
```

### Configuración inicial
```py
tile_id = (2, 2, "D") 
df = pd.read_csv(f"C:\\Users\\user\\Escritorio\\practicas\\tiles_B1-2\\tile_{tile_id[0]}_{tile_id[1]}_subtiles\\tile_{tile_id[0]}_{tile_id[1]}_{tile_id[2]}.csv") 

# Lista de genes de interés
genes_interes = ['Gfap','Nefm','Aldh1l1','Aqp4','Foxp3','Cd8a', 'Map2','Aldoc', 'Sparcl1',
                 'Flt3','Itgae','Itgax','Ptprc','Dct','Mitf','Pmel','Slc7a5','P2ry12','Cd3e','Cd4',
                 'Tmem119','Cd24a','Itgam','Ly6g','Adgre1','Nsg2','Rbfox3','Trem2','Klrg1',
                 'Nkx1-1','Mog','Olig1','Acta2','Fn1','Myh11','Pecam1']

genes_existentes = [g for g in genes_interes if g in df.columns]
if not genes_existentes:
    raise ValueError("Ninguno de los genes indicados está presente en el archivo.")

mask = df[genes_existentes].astype(int).sum(axis=1) > 0
df_filtrado = df[mask].copy()
```

### Función para agrupación
```py
def cluster_and_plot_scaled(df, genes_coloreados, distance_threshold=2.5, alpha_val=0.1):
        
    coords = df_filtrado[['x', 'y']].values
    gene_cols = df.drop(columns=['x', 'y']).select_dtypes(include=[np.number]).columns
    genes = df[gene_cols].values

    # Escalar por separado
    scaler_coords = StandardScaler()
    coords_scaled = scaler_coords.fit_transform(coords)

    scaler_genes = StandardScaler()
    genes_scaled = scaler_genes.fit_transform(genes)

    # Ponderar espacial
    peso_espacial = 15
    coords_scaled *= peso_espacial

    # Combinar
    features = np.hstack([coords_scaled, genes_scaled])

    # Clustering
    clustering = AgglomerativeClustering(
        n_clusters=None,
        distance_threshold=distance_threshold,
        linkage='average'
    )
    labels = clustering.fit_predict(features)
    df['cluster'] = labels

    print("Number of clusters:", len(set(labels)))

    # Plot
    fig, ax = plt.subplots(figsize=(8,8))
    ax.set_facecolor('white')

    # Genes coloreados
    for gene, color in genes_coloreados.items():
        if gene not in df.columns:
            continue
        sg = df[df[gene] == 1]
        ax.scatter(sg['x'], sg['y'], s=0.5, color=color, label=gene, alpha=0.5)

    # Contornos por cluster
    for cl in sorted(df['cluster'].unique()):
        cluster_points = df[df['cluster'] == cl]
        if len(cluster_points) < 8:
            continue
        
        # Gen dominante en el cluster
        gene_counts = {gene: cluster_points[gene].sum() for gene in gene_cols}
        gen_dominante = max(gene_counts, key=gene_counts.get)
        color_contorno = genes_coloreados.get(gen_dominante, 'black')

        # Contorno con alphashape
        try:
            puntos = cluster_points[['x', 'y']].values
            alpha_shape = alphashape.alphashape(puntos, alpha=alpha_val)
            if alpha_shape.is_empty:
                continue
            if isinstance(alpha_shape, Polygon):
                xs, ys = alpha_shape.exterior.xy
                ax.plot(xs, ys, color=color_contorno, linewidth=1.2)
            else:
                for geom in alpha_shape.geoms:
                    xs, ys = geom.exterior.xy
                    ax.plot(xs, ys, color=color_contorno, linewidth=1.2)
        except Exception as e:
            print(f"Cluster {cl} alphashape failed: {e}")

    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

    plt.legend(markerscale=5, fontsize=6)
    plt.savefig(f"cluster_scaled_tile.png", dpi=300, transparent=True, bbox_inches='tight', pad_inches=0)
    plt.show()

    return df
```

### Ejecutar la función
```py
clustered = cluster_and_plot_scaled( 
    df_filtrado,
    genes_coloreados=genes_coloreados, 
    distance_threshold=7,
    alpha_val=0.03
)
```

### Resultados


<p align="center">
  <img width="400" height="400" src="https://github.com/user-attachments/assets/c0c62d20-8662-4a15-bf82-ec7b159b28b8" />
</p>

<p align="center">
  <img width="350" height="350" src="https://github.com/user-attachments/assets/0be44cb1-53ae-4c99-9512-d8fd240db053" />
</p>

Se muestra la segmentación y debajo la imágen del DAPI con los nucleos correspondientes a esa SUBTILE.


## Agrupación por tipo celular

#### Módulos:
```py
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering
import alphashape
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
```

#### Función
```py
def cluster_and_plot_subset(
    df,
    genes_coloreados,
    genes_interes,
    distance_threshold=6,
    alpha_val=0.05,
    peso_espacial=8,
    expr_thr=0.0,          # umbral de “expresión positiva”
    k_at_least=1,          # al menos k genes_interes positivos
    min_pts_contour=10,     # mínimo de puntos para dibujar contorno
    point_size=0.5
):
    # --- 1) Validación y filtrado de columnas ---
    genes_validos = [g for g in genes_interes if g in df.columns]
    if not genes_validos:
        raise ValueError("Ninguno de los genes_interes está en el DataFrame.")

    # --- 2) Filtrar FILAS: solo coordenadas con >=k genes por encima de expr_thr ---
    presences = (df[genes_validos] > expr_thr)
    mask = presences.sum(axis=1) >= k_at_least
    df_sub = df.loc[mask, ['x','y'] + genes_validos].copy()
    if df_sub.empty:
        raise ValueError("Tras el filtrado no quedan puntos. Ajusta expr_thr o k_at_least.")

    # --- 3) Preparar matrices ---
    coords = df_sub[['x','y']].values
    genes = df_sub[genes_validos].values

    # Escalado
    coords_scaled = StandardScaler().fit_transform(coords) * peso_espacial
    genes_scaled  = StandardScaler().fit_transform(genes)

    features = np.hstack([coords_scaled, genes_scaled])

    # --- 4) Clustering jerárquico ---
    clustering = AgglomerativeClustering(
        n_clusters=None,
        distance_threshold=distance_threshold,
        linkage='average'
    )
    labels = clustering.fit_predict(features)
    df_sub['cluster'] = labels
    print(f"Nº clusters {genes_interes}: {df_sub['cluster'].nunique()}")

    # --- 5) Plot: SOLO el subset ---
    fig, ax = plt.subplots(figsize=(8,8))
    ax.set_facecolor('white')

    # Puntos por gen (solo en subset)
    for gene in genes_validos:
        color = genes_coloreados.get(gene, 'black')
        sg = df_sub[df_sub[gene] > expr_thr]
        if not sg.empty:
            ax.scatter(sg['x'], sg['y'], s=point_size, color=color, label=gene, alpha=0.5)

    # Contornos por cluster (dominante entre genes_validos)
    for cl in sorted(df_sub['cluster'].unique()):
        cluster_points = df_sub[df_sub['cluster'] == cl]
        if len(cluster_points) < min_pts_contour:
            continue

        gene_counts = {g: (cluster_points[g] > expr_thr).sum() for g in genes_validos}
        gen_dominante = max(gene_counts, key=gene_counts.get)
        color_contorno = genes_coloreados.get(gen_dominante, 'black')

        try:
            pts = cluster_points[['x','y']].values
            a = alphashape.alphashape(pts, alpha=alpha_val)
            if a.is_empty:
                continue
            if isinstance(a, Polygon):
                xs, ys = a.exterior.xy
                ax.plot(xs, ys, color=color_contorno, linewidth=1.2)
            else:
                for geom in a.geoms:
                    xs, ys = geom.exterior.xy
                    ax.plot(xs, ys, color=color_contorno, linewidth=1.2)
        except Exception as e:
            print(f"Cluster {cl} alphashape error: {e}")

    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    plt.legend(markerscale=5, fontsize=6)
    fname = f"clusters_subset_{'_'.join(genes_validos)}.png"
    plt.savefig(fname, dpi=300, transparent=True, bbox_inches='tight', pad_inches=0)
    plt.show()

    return df_sub
```

#### Ejecución de cada tipo celular:
```py
# T (helper/reg), B y NK
genes_lymph = ['Cd4','Foxp3','Cd8a','Cd3e','Cd19', 'Klrg1', 'Nkx1-1', 'Ptprc']
df_lymph = cluster_and_plot_subset(
    df_filtrado,
    genes_coloreados=genes_coloreados,
    genes_interes=genes_lymph,
    distance_threshold=7,
    alpha_val=0.05,
    expr_thr=0.0, # para datos binarios
    min_pts_contour=5,
    k_at_least=1
)
```

Resultados:
```py
Nº clusters ['Cd4', 'Foxp3', 'Cd8a', 'Cd3e', 'Cd19', 'Klrg1', 'Nkx1-1', 'Ptprc']: 12
```
<img width="640" height="636" alt="image" src="https://github.com/user-attachments/assets/4e8b5380-b286-41dc-870b-7fac84832116" />


```py
# Melanoma
genes_mel = ['Dct','Mitf','Pmel','Slc7a5']
df_mel = cluster_and_plot_subset(
    df_filtrado,
    genes_coloreados=genes_coloreados,
    genes_interes=genes_mel,
    distance_threshold=15,
    alpha_val=0.03,
    expr_thr=0.0,
    peso_espacial=40,
    min_pts_contour=40
)
```

Resultados:

```py
Nº clusters ['Dct', 'Mitf', 'Pmel', 'Slc7a5']: 62
```
<img width="640" height="636" alt="image" src="https://github.com/user-attachments/assets/fc94ad14-67bd-4c0b-bbe8-209abd835e39" />


```py
# Neuronas
genes_neu = ['Map2','Nefm','Nsg2','Rbfox3', 'Trem2']
df_neu = cluster_and_plot_subset(
    df_filtrado,
    genes_coloreados=genes_coloreados,
    genes_interes=genes_neu,
    distance_threshold=30,
    alpha_val=0.03,
    peso_espacial=55,
    min_pts_contour=10
)
```

Resultados:
```py
Nº clusters ['Map2', 'Nefm', 'Nsg2', 'Rbfox3', 'Trem2']: 31
```
<img width="640" height="636" alt="image" src="https://github.com/user-attachments/assets/a959abc8-8271-4005-bc43-490feaa6a014" />


```py
# Astrocitos
genes_ast = ['Gfap','Aldh1l1','Aldoc','Aqp4', 'Sparcl1']
df_ast = cluster_and_plot_subset(
    df_filtrado,
    genes_coloreados=genes_coloreados,
    genes_interes=genes_ast,
    distance_threshold=7,
    alpha_val=0.03,
    expr_thr=0.0,
    k_at_least=1
)
```

Resultados:
```py
Nº clusters ['Gfap', 'Aldh1l1', 'Aldoc', 'Aqp4', 'Sparcl1']: 39
```
<img width="640" height="636" alt="image" src="https://github.com/user-attachments/assets/4c571b9e-d1aa-4586-b218-c168c33fb3dc" />


```py
# Dendritic
genes_den = ['Flt3','Itgae','Itgax']
df_den = cluster_and_plot_subset(
    df_filtrado,
    genes_coloreados=genes_coloreados,
    genes_interes=genes_den,
    distance_threshold=8,
    alpha_val=0.01,
    peso_espacial=15,
    min_pts_contour=3
)
```

Resultados:
```py
Nº clusters ['Flt3', 'Itgae', 'Itgax']: 19
```
<img width="640" height="636" alt="image" src="https://github.com/user-attachments/assets/47f3163b-a12c-44c4-a741-67c7c02fe1b7" />


```py
# Microglia
genes_mic = ['P2ry12','Tmem119']
df_mic = cluster_and_plot_subset(
    df_filtrado,
    genes_coloreados=genes_coloreados,
    genes_interes=genes_mic,
    distance_threshold=7,
    alpha_val=0.01,
    peso_espacial=10,
    min_pts_contour=5
)
```

Resultados:
```py
Nº clusters ['P2ry12', 'Tmem119']: 12
```
<img width="640" height="636" alt="image" src="https://github.com/user-attachments/assets/18042ba8-cdb9-4b61-b4b5-65b03ccbf401" />


```py
# Mieloides
genes_miel = ['Cd24a','Itgam', 'Ly6g', 'Adgre1']
df_miel = cluster_and_plot_subset(
    df_filtrado,
    genes_coloreados=genes_coloreados,
    genes_interes=genes_miel,
    distance_threshold=8,
    alpha_val=0.02,
    peso_espacial=15,
    min_pts_contour=5
)
```

Resultados:
```py
Nº clusters ['Cd24a', 'Itgam', 'Ly6g', 'Adgre1']: 26
```
<img width="640" height="636" alt="image" src="https://github.com/user-attachments/assets/248458c5-c30e-40e8-9701-f5beea54bff7" />


```py
# Oligodendrocytes
genes_oli = ['Mog','Olig1']
df_oli = cluster_and_plot_subset(
    df_filtrado,
    genes_coloreados=genes_coloreados,
    genes_interes=genes_oli,
    distance_threshold=7,
    alpha_val=0.03,
    expr_thr=0.0,
    k_at_least=1
)
```
Resultados:
```py
Nº clusters ['Mog', 'Olig1']: 8
```
<img width="640" height="636" alt="image" src="https://github.com/user-attachments/assets/b7d97f92-6d31-4323-b022-d9a2ccd2b5c6" />


```py
# Vasculature
genes_vas = ['Acta2','Fn1', 'Myh11', 'Pecam1']
df_vas = cluster_and_plot_subset(
    df_filtrado,
    genes_coloreados=genes_coloreados,
    genes_interes=genes_vas,
    distance_threshold=6,
    alpha_val=0.02,
    peso_espacial=15,
    min_pts_contour=6
)
```

Resultados:
```py
Nº clusters ['Acta2', 'Fn1', 'Myh11', 'Pecam1']: 55
```
<img width="640" height="636" alt="image" src="https://github.com/user-attachments/assets/f8a38c29-ca03-4e44-98eb-e145526a5d8e" />




