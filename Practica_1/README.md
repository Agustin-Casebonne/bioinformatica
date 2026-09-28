# Simulador del dogma central de la biología molecular

Práctica 1 de Bioinformática (EII, ULPGC).

Simula en consola el flujo de la información genética:
ADN → ADN (replicación), ADN → ARNm (transcripción) y ARNm → proteína (traducción).

## Ejecución
python practica1.py

Requiere Python 3. No usa librerías externas.

## Entrada de la secuencia
1. Generada aleatoriamente.
2. Escrita por consola (solo A, T, C, G, hebra codificante).
3. Leída desde un fichero de texto o FASTA.

## Qué representa
- **Replicación**: helicasa, primasa (cebadores de ARN), ADN polimerasa III y I,
  ligasa, cadena líder y cadena rezagada con fragmentos de Okazaki.
- **Transcripción**: ARN polimerasa, hebra molde y complementariedad ADN → ARN.
- **Traducción**: ribosoma, lectura por codones desde AUG hasta el codón de parada.

## Simplificaciones
- Los fragmentos de Okazaki son de 6 bases y los cebadores de 3, por claridad.
  En la realidad son mucho más largos.
