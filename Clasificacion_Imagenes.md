# Analisis de Imagen Biomedica - Clasificación de imágenes

## Clasificación de imágenes con Bag of Features/Words en MATLAB. Comparativa resultados (Average Accuracy) de entrenamiento y validación

**Objetivo**: realizar una comparativa de resultados entre entrenamiento y validación midiendo el Average Accuracy del clasificador de categorías basado en Bag of Features.

``` splus
% Parametros
blockWidths = [32, 64, 96, 128];
vocabSizes = [5, 50, 250, 500];
results = [];   % guardar resultados

% Preparación de datos 
unzip('MerchData.zip');
imds = imageDatastore('MerchData','IncludeSubfolders',true,'LabelSource','foldernames');

% Separar train/validation sets
[trainingSet, validationSet] = splitEachLabel(imds, 0.6, 'randomize');

for bw = blockWidths
    for vs = vocabSizes
        fprintf('Training with BlockWidth=%d, VocabularySize=%d...\n', bw, vs);

        % Crear bag of features con los parametros
        bag = bagOfFeatures(trainingSet, 'VocabularySize', vs, 'BlockWidth', bw);

        % Train classifier
        categoryClassifier = trainImageCategoryClassifier(trainingSet, bag);

        % Evaluación entrenamiento y validación
        confTrain = evaluate(categoryClassifier, trainingSet);
        confVal   = evaluate(categoryClassifier, validationSet);

        accTrain = mean(diag(confTrain));
        accVal   = mean(diag(confVal));

        % Guardar resultados
        results = [results; {bw, vs, accTrain, accVal}];
    end
end

% Presentación de resultados
resultsTable = cell2table(results, ...
    'VariableNames', {'BlockWidth', 'VocabularySize', 'TrainAccuracy', 'ValidationAccuracy'});

disp(resultsTable)
```

**Resultados:**

| BlockWidth | VocabularySize | TrainAccuracy | ValidationAccuracy |
|------------|----------------|---------------|--------------------|
| 32         | 5              | 0.48889       | 0.5                |
| 32         | 50             | 0.8           | 0.7                |
| 32         | 250            | 0.97778       | 0.76667            |
| 32         | 500            | 1             | 0.9                |
| 64         | 5              | 0.53333       | 0.4                |
| 64         | 50             | 0.91111       | 0.83333            |
| 64         | 250            | 0.97778       | 0.93333            |
| 64         | 500            | 1             | 0.96667            |
| 96         | 5              | 0.6           | 0.53333            |
| 96         | 50             | 0.97778       | 0.93333            |
| 96         | 250            | 1             | 1                  |
| 96         | 500            | 1             | 0.96667            |
| 128        | 5              | 0.55556       | 0.43333            |
| 128        | 50             | 0.91111       | 0.86667            |
| 128        | 250            | 1             | 1                  |
| 128        | 500            | 1             | 1                  |

Graficamos los resulados:

``` splus
% Comparativa de precisión media en entrenamiento y validación
figure
tiledlayout(1,2)

% Subgráfica 1: Precisión en entrenamiento
nexttile
for i = 1:numel(blockWidths)
    idx = resultsTable.BlockWidth == blockWidths(i);
    plot(resultsTable.VocabularySize(idx), resultsTable.TrainAccuracy(idx), ...
        '-o', 'DisplayName', sprintf('BlockWidth=%d', blockWidths(i)));
    hold on
end
xlabel('Vocabulary Size')
ylabel('Train Accuracy')
title('Train Accuracy vs Vocabulary Size')
legend show

% Subgráfica 2: Precisión en validación
nexttile
for i = 1:numel(blockWidths)
    idx = resultsTable.BlockWidth == blockWidths(i);
    plot(resultsTable.VocabularySize(idx), resultsTable.ValidationAccuracy(idx), ...
        '-o', 'DisplayName', sprintf('BlockWidth=%d', blockWidths(i)));
    hold on
end
xlabel('Vocabulary Size')
ylabel('Validation Accuracy')
title('Validation Accuracy vs Vocabulary Size')
legend show
```

![](ex1.png)

<div style="page-break-before: always;"></div>

**Observaciones:**

> Se puede ver que `BlockWidth` y `VocabularySize` son parámetros clave, ya que la precisión varía bastante cuando los cambiamos. 
> 
> En términos generales, podemos observar que el tamaño de los bloques usado para extraer descriptores locales (`BlockWidth`) que mejor funciona para este dataset es el de **96** pixeles. Mientras que 32 es el que peor realiza la tarea tanto en entrenamiento como en validacion. Esto se debe a que al ser tan pequeño, está capturando detalles muy finos y locales, viendose afectado por el ruido y sobreentrenando, con lo que vemos una decaida en la validación.
>
> Por otro lado, con respecto al numero de palabras o clusters de descriptores locales que forman el vocabulario visual del modelo (`VocabularySize`), vemos que con **50** palabras todos los modelos mejoran drásticamente, especialmente el de tamaño 96 de bloque. A partir de ahí la mejora es ínfima. Para los tamaños de bloque grande se llega a un máximo en 250, que posteriormente disminuye. El tamaño 64, entra casi en plateau y el tamaño más pequeño continúa mejorando, muy levemente según se aumenta el número de palabras.


## Visualizar los 5 mejores “visual words” de cada categoría

**Objetivo**: El objetivo de este ejercicio es identificar y visualizar los 5 visual words más frecuentes en cada categoría de imágenes, utilizando Bag of Features con parámetros BlockWidth = 128 y VocabularySize = 500.

*Paso 1*. Preparación de datos, definición de parámetros y generación del Bag of Features:
``` splus
unzip('MerchData.zip');
imds = imageDatastore('MerchData', 'IncludeSubfolders', true, 'LabelSource', 'foldernames');

blockWidth  = 128;
vocabSize   = 500;

bag = bagOfFeatures(imds, 'VocabularySize', vocabSize, 'BlockWidth', blockWidth);

categories = unique(imds.Labels);
numCategories = numel(categories);
```
Esto genera el vocabulario visual necesario para codificar las imágenes.

*Paso 2*. Preparación de la carpeta de resultados
``` splus
outputDir = 'visual_words_output';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
```

*Paso 3*. Procesamiento por categoría

Para cada categoría:
- Se filtran las imágenes correspondientes.
- Se codifican todas las imágenes para obtener los índices de visual words presentes.
- Se identifican los 5 visual words más frecuentes.
- Se busca una imagen que contenga los 5 visual words seleccionados.
- Se simulan las posiciones de los visual words (por limitación de la versión de MATLAB).
- Se recortan y visualizan los parches asociados a cada visual word, guardando la figura resultante.

``` splus
for c = 1:numCategories
    % Filtrado, codificación y conteo de visual words
    category = categories(c);
    fprintf('\nProcesando categoría: %s\n', string(category));

    % Filtrar imágenes de la categoría
    idx = imds.Labels == category;
    idxFiles = find(idx);

    % Almacenar todos los índices de visual words
    allWords = [];

    % Codificar todas las imágenes de la categoría
    for i = 1:numel(idxFiles)
        I = readimage(imds, idxFiles(i));
        [~, words] = encode(bag, I);

        % Manejar distintos tipos de salida de 'words'
        if isa(words, 'vision.internal.visualWords')
            words = words.WordIndex;
        elseif iscell(words)
            words = cell2mat(words);
        end

        allWords = [allWords; double(words(:))];
    end

    % Calcular 5 visual words más frecuentes
    counts = histcounts(allWords, 1:vocabSize+1);
    [~, sortedIdx] = sort(counts, 'descend');
    top5 = sortedIdx(1:5);

    fprintf('Top 5 visual words para %s: %s\n', string(category), num2str(top5));

    % Selección de imagen que contenga los 5 visual words
    found = false;
    while ~found
        randIdx = randsample(idxFiles, 1);
        I = readimage(imds, randIdx);
        [~, words] = encode(bag, I);

        if isa(words, 'vision.internal.visualWords')
            words = words.WordIndex;
        elseif iscell(words)
            words = cell2mat(words);
        end
        words = double(words(:));

        if all(ismember(top5, words))
            found = true;
        end
    end

    % Simulación de posiciones y visualización de parches
    numWords = length(words);
    locations = rand(numWords, 2) .* [size(I,2), size(I,1)];

    figure('Name', sprintf('Categoría: %s', string(category)), 'NumberTitle', 'off');
    sgtitle(sprintf('%s: 5 Visual Words Más Comunes', string(category)));

    for j = 1:5
        wordIdx = top5(j);

        % Localizar una ocurrencia de ese visual word
        loc = locations(find(words == wordIdx, 1), :);
        cx = round(loc(1)); cy = round(loc(2));

        % Recorte alrededor del centro del patch
        halfBW = blockWidth / 2;
        x1 = max(cx - halfBW, 1);
        y1 = max(cy - halfBW, 1);
        x2 = min(cx + halfBW, size(I,2));
        y2 = min(cy + halfBW, size(I,1));

        patchImg = I(y1:y2, x1:x2, :);

        subplot(2, 3, j); % 2 filas × 3 columnas
        imshow(patchImg);
        title(sprintf('Patch %d (Word %d)', j, wordIdx), 'FontSize', 8);
    end

    % Guardar imagen
    saveas(gcf, fullfile(outputDir, sprintf('%s_top5.png', string(category))));
end

fprintf('\n Las imágenes se guardaron en la carpeta "%s"\n', outputDir);
```

<div style="page-break-before: always;"></div>

**Resultado:**

![](visual_words_output/MathWorks%20Cap_top5.png)
![](visual_words_output/MathWorks%20Cube_top5.png)
![](visual_words_output/MathWorks%20Playing%20Cards_top5.png)
![](visual_words_output/MathWorks%20Screwdriver_top5.png)
![](visual_words_output/MathWorks%20Torch_top5.png)

**Observaciones:**

> Al mostrar las 5 palabras más comunes de las 5 fotos vemos que algunas son sólo fondo. Esto se debe seguramente a que hemos fijado el tamaño de los boques y, como algunos objetos son más grandes que otros, el tamaño de estos bloques puede no ser el más eficiente para todos los casos. En las imágenes de la gorra y el cubo, las 5 palabras son parte del objeto; pero, en las imágenes de la baraja, el destornillador y la linterna, hay al menos 2 palabras que son solo fondo. Esto es un problema porque el modelo está relacionando ese color de fondo con el objeto (todos los fondos son distintos tonos de azul), lo que puede producir sobre-aprendizaje.
Lo mejor sería adaptar el tamaño del bloque a cada imagen en particular.