## CytoMAP

### Install
Installing the standalone [version](https://gitlab.com/gernerlab/cytomap/-/blob/master/StandaloneInstaller/CytoMAP_Installer_WindowsV1.4.1.exe)

### Código para transformar el df segmentado a formato CytoMAP
Para usar la segmentación de `AgglomerativeClustering` en CytoMAP es necesario exportar un CSV donde cada fila sea una célula con: 

- CellID → identificador único por célula 
- X_position, Y_position → coordenadas en µm (CytoMAP las espera en micras, así que si tus coordenadas están en píxeles, hay que convertir) 
- ClusterID → número de cluster que sale del agglomerative clustering 
- CellType → por ejemplo, el gen dominante en esa célula 
- Marcadores → opcional, intensidades o 0/1 para cada gen (para análisis adicionales en CytoMAP) 

:construction: En construcción... :construction:

```py

```

```py

```



```py

```


```py

```


```py

```


```py

```


```py

```
