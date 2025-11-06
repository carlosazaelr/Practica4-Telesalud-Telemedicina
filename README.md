# Reporte de Práctica 3: Introducción al manejo de señales biomédicas: Señal de ECG
Telesalud y Telemedicina  
Ingeniería Biomédica  
Departamento de Ingeniería Eléctrica y Electrónica, Tecnológico Nacional de México/IT Tijuana, Blvd. Alberto Limón Padilla s/n, Tijuana, C.P. 22454, B.C., México.  
Prof. Fortunato Ramírez Arzate  
Unidad 1: Introducción a la Telesalud y a la Telemedicina  
Ramirez Rodriguez Carlos Azael  
22212267  

## Introducción

El estudio de las señales bioeléctricas constituye una de las bases fundamentales en el análisis de sistemas biomédicos, ya que permite comprender el comportamiento fisiológico del cuerpo humano a través de la actividad eléctrica de diferentes órganos. Entre estas señales, la señal electrocardiográfica (ECG) es una de las más relevantes, pues refleja la actividad eléctrica del corazón y proporciona información esencial para el diagnóstico y monitoreo de patologías cardiovasculares.

En esta práctica se trabajará con el procesamiento digital de señales de ECG empleando **MATLAB** como herramienta de análisis. El enfoque principal estará en tres etapas fundamentales del tratamiento de señales:

- **Normalización**: permite ajustar la amplitud de la señal a un rango predefinido, lo cual facilita su comparación, almacenamiento y posterior procesamiento.
- **Amplificación**: mejora la visibilidad de la señal, resaltando componentes de interés que de otra forma serían difíciles de analizar debido a su baja magnitud.
- **Cuantificación**: consiste en representar la señal continua en niveles discretos de amplitud, emulando el proceso realizado en los sistemas de adquisición mediante convertidores analógico-digitales (ADC).

## Descarga/Adquisición de datos

Para iniciar con el desarrollo de la práctica, se obtuvo un archivo de una señal de ECG desde una base de datos en línea. La base de datos en cuestión fue la siguiente:

* **Señal de electrocardiograma (ECG)**. Utilizar una señal de la base de datos que se desee. Por ejemplo:

  * [PhysioNet Database](https://physionet.org/about/database/)  

El archivo utilizado para el desarrollo de la práctica fue el siguiente: [ecgdata.mat](p_3/ecgdata.mat)

## Configuración y despliegue en Matlab

Una vez descargado nuestro archivo, se guardó en la carpeta de trabajo para abrirla desde la interfaz de **MATLAB**, tal como se muestra en la siguiente figura.  

<p align="center">
  <img src="imagenes/intro.png" alt="Espacio de trabajo" width="400">  
</p>
<p align="center"><em>Figura 1. Espacio de trabajo</em></p>  

Para realizar el despliegue de la figura de la señal ECG original, se utilizaron las siguientes líneas de código:

```matlab
% Parámetros de la señal
Fs = 2000;
Obits = 16;
ini = 170; % [s]
fin = 172; % [s]

% Despliegue de la figura
figure;
load ecgdata.mat
x = ecgdata;
x = x(Fs*ini:(Fs*fin)-1);
t = ini:1/Fs:fin-1/Fs;
plot(t, x);
```
Donde: Fs = Frecuencia de muestreo, Obits = # de bits, ini = Segundo de inicio, fin = Segundo final.

Una vez ejecutado el código, se obtuvo la señal representada en la Figura 2.  

<p align="center">
  <img src="imagenes/original.jpg" alt="ECG original" width="800">  
</p>
<p align="center"><em>Figura 2. Señal ECG original </em></p>

En orden para preparar la señal para su posterior cuantificación, se necesita trasladar, normalizar y amplificar. Para ello, se utilizó la siguiente sección de código.

```matlab
% Trasladar
x_t = x + abs(min(x));

% Normalizar
x_n = x_t / max(x_t);

% Amplificar
x_a = x_n * 2^Obits-1;
```
Una vez ejecutado, se expresó este proceso en una serie de gráficos, la cual se puede observar en la Figura 3.

<p align="center">
  <img src="imagenes/fase2.jpg" alt="ECG amplificada" width="800">  
</p>
<p align="center"><em>Figura 3. Señal ECG trasladada, normalizada y amplificada </em></p>

## Cuantificación de la señal

Para este paso, se definió una función de nombre **cuantificador**, la cual tendrá como entradas `Xin` (Señal de entrada en forma de vector), `Nbits` (Número de bits del cuantificador) y `Obtis` (Número de bits de la señal). La función correspondiente se definió de la siguiente forma:

```matlab
function Yout = cuantificador_3(Xin, Nbits, Obits)
    P = (2^Obits/(2^Nbits))-1 : 2^Obits/(2^Nbits) : 2^Obits-1;
    Yout = quantiz (Xin,P);
end
```

## Generación y análisis de fragmentos

Para esta sección del desarrollo, se siguieron los siguientes puntos:
1. Elegir un fragmento de la señal de entre **2 y 4 segundos** para el análisis.
2. Generar varias versiones de la señal de ECG original usando diferentes valores de `Nbits`.
3. Desplegar de manera simultánea la señal original y las señales cuantificadas (utilizar el comando `subplot`) y comparar los resultados.
4. Explicar cómo afectan los diferentes valores del parámetro `Nbits`.

Para ello, se utilizaron las siguientes instrucciones:

```matlab
figure;
for n = 1: 1 :Obits
    subplot(4,4,n);
    x_c = cuantificador_3(x_a,n,Obits);
    plot(t,x_c)
    title(strcat(num2str(n),'-bit'))
    xlabel('Tiempo');
    ylabel('Amplitud');
    xlim([ini fin]);
end
```

Al ejecutar el ciclo del listado anterior, se obtiene un mosaico como el de la Figura 4.

<p align="center">
  <img src="imagenes/fase3.jpg" alt="ECG mosaico" width="800">  
</p>
<p align="center"><em>Figura 4. Señal ECG cuantificada a diferentes valores de Nbits </em></p>

## Análisis de resultados

Para realizar el análisis de resultados, se acomodaron los datos en la siguiente tabla, junto a las observaciones pertinentes.

| Profundidad de bit ó Nbits (niveles de cuantificación) | Calidad (MB/B/R/M/MM) | Observaciones                    |
| ------------------------------------------------------ | --------------------- | -------------------------------- |
| 1 (2)                                                  | MM                    | La señal toma únicamene valores de 0 o 1, lo cual no nos deja apreciar bien su forma, pero funciona como un detector de picos. |
| 2 (4)                                                  | MM                    | Igual que en el nivel anterior, pero se logra apreciar la forma de los picos. |
| 3 (8)                                                  | MM                    | Al igual que el nivel anterior, con la diferencia de que apenas se alcanzan a notar los niveles inferiores de la señal.                                |
| 4 (16)                                                 | MM                    | Se puede empezar a distinguir la forma de la señal pero con una muy baja calidad.                                |
| 5 (32)                                                 | M                     | Un poco mejor que el nivel anterior, aunque aún se nota la pérdida de datos.                                |
| 6 (64)                                                 | M                     | Igual que en el caso anterior con una cierta mejora.                                |
| 7 (128)                                                | R                     | El mejor nivel en cuanto a calidad y peso en bits.                                |
| 8 (256)                                                | R                     | Cierta mejora con respecto al anterior.                                |
| 9 (512)                                                | B                     | Se puede apreciar una imagen más limpia.                                |
| 10 (1024)                                              | B                     | Poco cambio con respecto a la anterior                                |
| 11 (2048)                                              | MB                    | Imagen muy similar a la original.                                |
| 12 (4096)                                              | MB                    | Sin cambios significativos.                                |
| 13 (8192)                                              | MB                    | Sin cambios significativos.                                |
| 14 (16384)                                             | MB                    | Sin cambios significativos.                                |
| 15 (32768)                                             | MB                    | Sin cambios significativos.                                |
| 16 (65536)                                             | MB                    | La mejor calidad. |

**¿Cuál es la profundidad de bit adecuada para la señal ECG utilizada en esta práctica?**  
En mi opinión, el valor de 7 en escala de profundidad de bits. En este nivel se captura de una buena manera las formas que se deberían poder identificar hablando de una señal de ECG, siendo un peso en bits bastante menor al de la imagen original. Sin embargo, la señal de profundidad de bit de 1 puede ser útil si solo se desean contabilizar los picos de la señal.

## Conclusión

En esta práctica se logró familiarizarse con el manejo de señales bioeléctricas de ECG utilizando MATLAB, aplicando etapas fundamentales de procesamiento digital como la **traslación, normalización, amplificación y cuantificación**. Se evidenció cómo cada una de estas etapas afecta la representación y análisis de la señal:

- La **traslación** permitió que todas las muestras de la señal fueran positivas, facilitando su procesamiento posterior.  
- La **normalización** ajustó la señal a un rango común, haciendo más fácil la comparación entre diferentes fragmentos y señales.  
- La **amplificación** incrementó la magnitud de la señal hasta el rango completo definido por el número de bits de la adquisición, resaltando detalles importantes.  
- La **cuantificación**, mediante diferentes niveles de Nbits, mostró cómo la precisión digital afecta la fidelidad de la señal: valores bajos resultan en una representación muy simplificada, mientras que valores mayores permiten conservar la forma de la señal con buena resolución.  

Se concluye que, para la señal ECG utilizada en esta práctica, una **profundidad de 7 bits** proporciona un balance adecuado entre calidad de la señal y eficiencia de almacenamiento, capturando los picos y las características esenciales de la señal sin generar un exceso de datos. Además, se observó que profundidades de bit menores pueden ser útiles para aplicaciones específicas como la detección de picos, mientras que profundidades mayores no ofrecen mejoras significativas perceptibles en la señal.  

En general, esta práctica permitió comprender la importancia de cada etapa del procesamiento de señales bioeléctricas y cómo las decisiones de cuantificación impactan directamente en el análisis y la interpretación de datos biomédicos.

