# Plan de desarrollo: ARMS-PCR

## Propósito

Incorporar al paquete `rflpSNP` el diseño *in silico* de ensayos
tetra-primer ARMS-PCR para genotipificar SNPs cuando PCR-RFLP no es viable
porque el cambio alélico no crea ni elimina un sitio de restricción útil.

La primera versión se centrará en SNPs bialélicos y en una única reacción con
cuatro cebadores. El resultado esperado será interpretable mediante
electroforesis en gel de agarosa para los tres genotipos: homocigoto de
referencia, heterocigoto y homocigoto alternativo.

## Alcance y límites

Incluido en esta iniciativa:

- Lectura de una secuencia FASTA y localización del SNP con la infraestructura
  actual del paquete.
- Diseño, evaluación y ranking de conjuntos de cuatro cebadores para
  tetra-primer ARMS-PCR.
- Exportación de un reporte detallado en texto.
- Simulación *in silico* de los productos ARMS-PCR y de los patrones de gel.
- Documentación, ejemplos y pruebas automatizadas.

Excluido explícitamente:

- Modificar, importar, adaptar o revisar código de las ramas PIRA-PCR.
- PIRA-PCR, PCR-RFLP con desajuste inducido, multiplexado de varios SNPs y
  variantes multialélicas.
- Validación experimental, predicción de especificidad genómica mediante BLAST
  o uso diagnóstico clínico.

Las funciones ARMS-PCR deben ser independientes de las rutas PIRA-PCR. Solo
podrán reutilizarse utilidades ya presentes en la rama principal que sean
generales, estables y pertinentes, como cálculo de Tm, contenido GC y cribado
heurístico de dímeros/hairpins.

## Decisión técnica: tetra-primer ARMS-PCR

Cada ensayo contendrá cuatro cebadores:

1. Dos **externos**, que flanquean el SNP y producen un amplicón de control.
2. Un cebador **interno específico del alelo de referencia**.
3. Un cebador **interno específico del alelo alternativo**.

Los dos internos terminan en el SNP en su extremo 3'. Cada uno incluirá una
segunda discrepancia deliberada, inicialmente configurable en la posición -2
o -3 desde 3', para reducir la extensión sobre el alelo no objetivo. La
selección de esa discrepancia se considerará una recomendación de diseño y no
una garantía experimental.

Se favorecerán diseños en que:

- ambos amplicones alélicos sean distinguibles entre sí;
- ambos sean distinguibles del amplicón externo de control;
- los cuatro cebadores tengan propiedades fisicoquímicas compatibles;
- el conjunto completo tenga bajo riesgo de dímeros y hairpins;
- cada cebador interno sea refractario, según la regla de discrepancias 3', al
  alelo no objetivo.

## Flujo de usuario objetivo

```text
FASTA de referencia
        |
        v
Secuencia de flanqueo
        |
        v
Coordenada y alelo de referencia del SNP
        |
        +-- alelo alternativo indicado explícitamente por el usuario
        |
        v
Diseño y ranking de sets ARMS (4 cebadores)
        |
        +--> reporte TXT
        |
        v
Simulación de PCR ARMS para los 3 genotipos
        |
        v
Simulación de gel de agarosa
```

El alelo alternativo será siempre un argumento explícito de las funciones
ARMS-PCR de la primera versión. No se inferirá solamente a partir de la
secuencia flanqueante.

## Interfaz propuesta

Los nombres son parte del contrato inicial y podrán ajustarse antes de la
implementación, sin cambiar el alcance funcional.

### `design_arms_primers()`

Entradas principales:

- `gene_seq`: `DNAString` con la secuencia de referencia.
- `snp_pos`: coordenada 1-based del SNP.
- `ref_allele` y `alt_allele`: bases de los dos alelos.
- Restricciones de Tm, GC, longitudes, tamaños de amplicón, separación mínima
  de bandas y parámetros de la segunda discrepancia.

Salida: objeto S3 `rflp_arms_primers` con `top`, `best_set`, número total de
sets válidos, parámetros usados y diagnósticos de filtrado. `top` devolverá
por defecto los cinco mejores sets.

Cada fila describirá los cuatro cebadores, sus secuencias 5'→3', coordenadas,
orientaciones, Tm, GC, cambios deliberados, tamaños esperados de los tres
productos y una puntuación reproducible.

### `export_arms_primers_txt()`

Exportará el set recomendado, los cinco mejores candidatos, los parámetros
empleados, las bandas predichas y advertencias de validación experimental.

### `simulate_arms_pcr()`

Recibirá el set de cebadores y generará los productos esperados para
`ref/ref`, `ref/alt` y `alt/alt`. Debe verificar explícitamente la orientación,
las coordenadas y la especificidad de los extremos 3' de los internos; no
debe reutilizar sin cambios el límite genérico de discrepancias de
`simulate_pcr()`.

### `simulate_arms_gel()`

Construirá un `ggplot` con marcador, control externo y bandas diagnósticas,
mostrando una pista para cada genotipo esperado. No modificará el
comportamiento actual de `simulate_gel()`.

### `run_arms_pcr_pipeline()`

Encadenará `locate_snp()`, `design_arms_primers()`, la exportación opcional,
`simulate_arms_pcr()` y `simulate_arms_gel()`. Devolverá los resultados de
cada etapa para permitir inspección y reproducibilidad.

## Etapas de trabajo

### Etapa 0 — Base verificable

- Corregir el metadato de `DESCRIPTION` que impide que `R CMD check` complete
  la comprobación actual.
- Crear la estructura de pruebas y los casos sintéticos independientes de
  datos externos.
- Confirmar las convenciones existentes de coordenadas, orientación y objetos
  `DNAString` que usarán las nuevas funciones.

Criterio de salida: las comprobaciones del paquete pasan y hay pruebas que
describen los resultados esperados antes de implementar el algoritmo.

### Etapa 1 — Contrato y generador de candidatos ARMS

- Definir el formato exacto de `rflp_arms_primers`.
- Generar candidatos externos e internos en ambas orientaciones posibles.
- Forzar el SNP como nucleótido terminal 3' de cada interno.
- Enumerar y registrar la segunda discrepancia candidata en -2/-3.

Criterio de salida: se generan conjuntos geométricamente válidos de cuatro
cebadores, sin hacer todavía una recomendación final.

### Etapa 2 — Filtrado, especificidad y ranking

- Aplicar límites de longitud, Tm, GC, hairpins y dímeros individuales y
  cruzados para los cuatro cebadores.
- Validar la complementariedad perfecta contra el alelo objetivo y la doble
  discrepancia terminal contra el no objetivo.
- Filtrar por tamaños y separación de los tres amplicones.
- Puntuar y devolver los cinco mejores sets, con diagnósticos de exclusión.

Criterio de salida: `design_arms_primers()` produce resultados reproducibles y
explica por qué un SNP no obtuvo un set válido.

### Etapa 3 — Reporte y simulación

- Implementar `export_arms_primers_txt()`.
- Implementar `simulate_arms_pcr()` para los tres genotipos.
- Implementar `simulate_arms_gel()` y mapas de bandas interpretables.
- Añadir la función orquestadora `run_arms_pcr_pipeline()`.

Criterio de salida: un usuario puede ir de FASTA y flanqueo a un reporte y un
gel virtual con una sola llamada, conservando las salidas intermedias.

### Etapa 4 — Calidad y documentación

- Añadir documentación roxygen y páginas de ayuda.
- Actualizar `README.md` con un ejemplo completo y limitaciones.
- Añadir pruebas de regresión y ejecutar `R CMD check`.
- Revisar compatibilidad de los objetos nuevos con las versiones de R
soportadas por la integración continua.

Criterio de salida: la característica queda documentada, probada y lista para
revisión en una pull request desde `arms_pcr`.

## Matriz mínima de pruebas

- Un SNP bialélico con diseño válido en cada orientación.
- Transiciones y transversiones.
- Genotipos `ref/ref`, `ref/alt` y `alt/alt` con bandas esperadas inequívocas.
- Regiones cercanas a los extremos de la secuencia.
- Candidatos descartados por Tm, GC, dímeros, hairpins y bandas solapadas.
- Casos donde los internos amplifican indebidamente el alelo no objetivo.
- Conservación del comportamiento de todas las funciones PCR-RFLP actuales.

## Límites de interpretación experimental

El módulo ofrecerá diseños candidatos, no validación experimental. Antes de
sintetizar los oligonucleótidos se deberán confirmar especificidad genómica,
estructuras secundarias con una herramienta termodinámica especializada y
condiciones de PCR. ARMS-PCR requiere una polimerasa sin actividad correctora
3'→5'; la elección final de la segunda discrepancia y la temperatura de
alineamiento deben validarse en el laboratorio.

## Continuidad entre tareas

Este archivo es la fuente de verdad para el desarrollo de ARMS-PCR. Cada nueva
tarea deberá indicar la etapa a continuar, la rama `arms_pcr` y el requisito
de no tocar ni consultar PIRA-PCR. Los cambios se harán en etapas pequeñas,
con pruebas y una revisión del estado del repositorio antes de iniciar la
siguiente.
