DEL ADN a las proteínas

Autores:

Jaime Rivero Santana 
Agustín Darío Casebonne


Repositorio: 

https://github.com/Agustin-Casebonne/bioinformatica/tree/main/Practica%20Del%20ADN%20a%20la%20proteína 

1. Replicación del ADN
La replicación del ADN es un proceso semiconservativo en el que cada molécula de ADN resultante conserva una hebra de la molécula original y sintetiza una nueva hebra complementaria.
Partiendo de la secuencia:
5' - ATG CCG TTA GCT - 3'
3' - TAC GGC AAT CGA - 5'

    Tras una ronda de replicación se obtienen dos moléculas de ADN, cada una formada por una hebra original y una nueva:
    Molécula 1: 5' - ATG CCG TTA GCT - 3' (original)
                        3' - TAC GGC AAT CGA - 5' (nueva)

    Molécula 2: 5' - ATG CCG TTA GCT - 3' (nueva)
                            3' - TAC GGC AAT CGA - 5' (original)

    Durante el proceso, la helicasa separa las dos hebras, la primasa coloca los cebadores, la ADN polimerasa sintetiza las nuevas cadenas y la ligasa une los fragmentos de la hebra retardada.
    Mediante Biopython se obtuvo automáticamente la hebra complementaria de la secuencia ATGCCGTTAGCT, obteniéndose TACGGCAATCGA, coincidiendo con el resultado calculado manualmente. Para ello se utilizó Seq, que permite trabajar con secuencias biológicas, y la función complement(), que genera la secuencia complementaria.
    Si un error del ADN polimerasa no es corregido, puede convertirse en una mutación. Dependiendo de la región afectada, esta mutación podría modificar posteriormente el ARN o la proteína producida y alterar su función.

2. Transcripción del ADN a ARN
La transcripción consiste en utilizar una de las hebras de ADN como molde para sintetizar una molécula de ARN complementaria en dirección 5' → 3'. En este ejercicio, la hebra molde es:
3' - TAC GGA CTT ACG - 5'
A partir de ella se obtiene el siguiente ARNm:
5' - AUG CCU GAA UGC - 3'
La ARN polimerasa utiliza la hebra molde como referencia para incorporar los nucleótidos complementarios, utilizando uracilo (U) en lugar de timina (T). El promotor es una zona del ADN que indica al ARN polimerasa dónde comienza la transcripción, mientras que la región codificante contiene la información que se utilizará para producir la proteína. En este ejercicio no se muestra el promotor, sólo una parte de la secuencia de ADN.
Como extensión con Biopython, la secuencia de ADN se almacenó en un archivo FASTA y se leyó mediante SeqIO.read(). Posteriormente, complement() permitió obtener la secuencia complementaria y transcribe() convertirla en ARN, obteniéndose AUGCCUGAAUGC, que coincide con el resultado calculado manualmente.

3. Traducción del ARNm a proteína
La traducción consiste en convertir la información del ARNm en una secuencia de aminoácidos. La secuencia utilizada fue:
5' - AUG UAU GCU UAA - 3'
El codón AUG indica el inicio de la traducción, UAU codifica para tirosina y GCU para alanina. Finalmente, UAA indica que la traducción debe terminar. Por tanto, la proteína obtenida es:
Metionina - Tirosina - Alanina
Mediante Biopython, la función translate() produjo MYA*, donde M, Y y A representan los tres aminoácidos y * representa el codón de parada.
Si AUG se sustituye por GUG, se pierde el codón de inicio habitual y la traducción puede comenzar de manera incorrecta. Si desaparece el codón UAA, el ribosoma continuará leyendo la secuencia hasta encontrar otro codón de parada.

4. Splicing alternativo
El splicing alternativo permite que un mismo gen genere diferentes moléculas de ARNm mediante la selección de distintas combinaciones de exones. Considerando un gen formado por cinco exones (1-2-3-4-5), se pueden obtener, por ejemplo, los siguientes transcritos:
Transcrito 1: 1-2-4-5
Transcrito 2: 1-3-5
Al utilizar diferentes combinaciones de exones, las secuencias de los ARNm son diferentes y pueden dar lugar a proteínas con distintas secuencias, longitudes o funciones. En el programa desarrollado, los dos transcritos produjeron las proteínas MKPGKFG y MKFPG, respectivamente.
Este mecanismo aumenta la diversidad proteica porque un mismo gen puede producir diferentes transcritos y, por tanto, diferentes proteínas, sin necesidad de disponer de un gen independiente para cada una.
Como ejemplo real se utilizó el gen humano FGFR2, que presenta múltiples transcritos en Ensembl. Estos transcritos pueden presentar diferencias en los exones utilizados y en la longitud de las proteínas producidas. Estas variaciones pueden modificar determinadas regiones de la proteína y, como consecuencia, sus propiedades o función.
Para la implementación se utilizó Seq de Biopython para representar las secuencias de los exones y translate() para traducir los diferentes transcritos a proteínas y observar las diferencias producidas por el splicing alternativo.

5. Introducción a las proteínas
La secuencia de aminoácidos considerada es:
Met – Ile – Ser – Gly – Val – Lys – His
El extremo N-terminal corresponde a la metionina (Met), que es el primer aminoácido de la cadena, mientras que el extremo C-terminal corresponde a la histidina (His), que es el último.
El orden de los aminoácidos influye directamente en la estructura final de la proteína, ya que determina las interacciones entre sus aminoácidos y, por tanto, cómo se pliega la cadena en una estructura tridimensional. Si una mutación cambiara un aminoácido hidrofóbico por uno hidrofílico en una región interna de la proteína, podrían modificarse las interacciones en esa zona y alterar el plegamiento. Esto podría afectar a la estructura y, en determinados casos, a la función de la proteína.
Como extensión bioinformática, se consultó el Protein Data Bank (PDB) para observar la estructura tridimensional de una proteína conocida. Se utilizó la estructura 1HHO, correspondiente a la hemoglobina humana, en la que se pueden observar principalmente estructuras secundarias en forma de hélices alfa. Esta observación permite relacionar la secuencia de aminoácidos con el plegamiento de la proteína. Una mutación puntual podría alterar determinadas interacciones entre aminoácidos y modificar la estructura tridimensional si afecta a una región importante para el plegamiento.

6. Actividad integradora: del ADN a la proteína
Para realizar esta actividad se utilizó una secuencia de ADN obtenida de la base de datos pública NCBI, concretamente del registro Homo sapiens insulin (INS), transcript variant 1, mRNA, con número de acceso NM_000207.3.
 En primer lugar, se descargó el registro completo en formato FASTA. Como el registro contiene el ARNm completo, se consultó su anotación y se identificó la región codificante (CDS), correspondiente a las posiciones 60..392. Para obtener únicamente esta región, se desarrolló el programa extraer_cds.py, que lee el archivo FASTA mediante SeqIO.read(), extrae las posiciones correspondientes a la CDS y genera un nuevo archivo adn.fasta con los 333 nucleótidos de la secuencia codificante.
En primer lugar, se realizó la replicación obteniendo la hebra complementaria. De esta forma se obtienen dos moléculas de ADN, cada una formada por una hebra original y una nueva.
A continuación, se utilizó la hebra complementaria como molde para realizar la transcripción.
Finalmente, el ARNm se tradujo utilizando el código genético, obteniéndose:
MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN*
El símbolo * representa un codón de parada. El programa pipeline.py integra los tres procesos y muestra al usuario las secuencias obtenidas en cada etapa. Para ello se utilizaron las funciones complement() para obtener la secuencia complementaria, transcribe() para obtener el ARN y translate() para traducir el ARNm a proteína.
Un error durante cualquiera de estos procesos puede afectar al resultado final. Una mutación en el ADN puede modificar un codón del ARNm y, en consecuencia, cambiar un aminoácido, introducir un codón de parada o alterar la proteína producida. Por ello, los errores que modifican la secuencia codificante pueden tener consecuencias sobre la estructura y función de la proteína.
