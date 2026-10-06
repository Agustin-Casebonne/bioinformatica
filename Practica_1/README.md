# Simulador del dogma central de la biología molecular

**Práctica 1 de Bioinformática**

**Escuela de Ingeniería Informática - Universidad de Las Palmas de Gran Canaria (ULPGC)**

## Descripción 

Este proyecto desarrolla un simulador educativo del dogma central de la biología molecular, cuyo objetivo es representar mediante un programa el flujo de la información genética a través de los principales procesos de la biología molecular:

ADN → ADN (replicación), ADN → ARNm (transcripción) y ARNm → proteína (traducción).

El programa está desarrollado en Python y funciona mediante una interfaz de consola. A partir de una secuencia de ADN, el simulador muestra de forma secuencial los principales acontecimientos de la replicación, la transcripción y la traducción.

La finalidad del programa es facilitar la comprensión de estos procesos mediante una representación computacional sencilla y visual de las moléculas, enzimas y reglas de complementariedad implicadas.

## Ejecución

El programa se ejecuta desde la carpeta que contiene el archivo practica1.py:
```bash
python practica1.py
```
## Requisitos 

- Python 3.
- No se requieren librerías externas. Se utiliza únicamente la biblioteca estándar de Python.
  
## Entrada de la secuencia

1. Generación aleatoria: el programa genera automáticamente una secuencia de ADN.

2. Entrada manual: el usuario introduce por consola la hebra codificante utilizando únicamente las bases A, T, C y G.
   
3. Lectura desde fichero: se puede proporcionar una secuencia mediante un fichero de texto o FASTA.

## Qué representa

- **Replicación**: helicasa, primasa (cebadores de ARN), ADN polimerasa III y I,
  ligasa, cadena líder y cadena rezagada con fragmentos de Okazaki.
- **Transcripción**: ARN polimerasa, hebra molde y complementariedad ADN → ARN.
- **Traducción**: ribosoma, lectura ARNm mediante codones desde el codón de inicio AUG hasta un codón de parada.

## Simplificaciones

El programa es una representación educativa de los procesos biológicos y utiliza algunas simplificaciones para facilitar su visualización:

- Los fragmentos de Okazaki se representan con una longitud de 6 bases.
- Los cebadores de ARN se representan con una longitud de 3 bases.

## Repositorio 

Enlace al repositorio: [Práctica 1](https://github.com/Agustin-Casebonne/bioinformatica/tree/main/Practica_1)

## Autores

Jaime Rivero Santana

Agustín Darío Casebonne
