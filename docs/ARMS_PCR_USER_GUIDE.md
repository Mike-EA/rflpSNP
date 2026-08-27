# Guía de usuario: tetra-primer ARMS-PCR

Esta guía describe el flujo ARMS-PCR de `rflpSNP` para SNPs bialélicos que no
cuentan con un cambio de sitio de restricción aprovechable. El resultado es un
conjunto de cuatro primers, productos esperados para los tres genotipos, un
reporte de texto, un gel virtual y mapas nucleotídicos reproducibles.

## Alcance y límites

El módulo diseña candidatos *in silico*. No confirma especificidad genómica,
no sustituye validación termodinámica especializada y no garantiza rendimiento
experimental. Use una polimerasa sin actividad correctora 3'→5' y valide en
laboratorio la temperatura de alineamiento y la segunda discrepancia.

## Flujo básico

```r
devtools::load_all(".")

gene_seq <- read_gene_fasta("referencia.fasta")
snp <- locate_snp(gene_seq, flank_seq = "SECUENCIA_5_PRIMA_ANTES_DEL_SNP")

design <- design_arms_primers(
  gene_seq,
  snp_pos = snp$snp_pos,
  ref_allele = snp$snp_base,
  alt_allele = "A"
)

design$best_set
```

`ref_allele` siempre debe coincidir con la base de la referencia en
`snp_pos`; `alt_allele` debe ser la otra base en la misma orientación de la
referencia FASTA.

## Interpretar el diseño

Un resultado válido incluye:

- dos primers externos que producen la banda de control;
- un interno específico del alelo de referencia;
- un interno específico del alelo alternativo;
- tres tamaños: control, referencia y alternativa.

`design$diagnostics` resume el filtrado individual por Tm, GC, auto-dímero y
hairpin. `design$exclusion_diagnostics` resume por qué los sets completos se
descartaron. Si `n_valid_sets` es cero, revise primero los pools internos:

```r
design$diagnostics
design$exclusion_diagnostics
```

## Exploración de regiones difíciles

Para una exploración rápida, `max_raw_candidates_per_pool` limita cuántos
candidatos se evalúan antes del cálculo termodinámico. Es útil en ventanas
grandes, pero no reemplaza una búsqueda exhaustiva para una recomendación
final.

El caso de regresión FTO rs8050136 incluido en el paquete usa:

```r
design <- design_arms_primers(
  gene_seq, snp_pos, ref_allele = "C", alt_allele = "A",
  outer_flank = 85L, outer_length_min = 20L, outer_length_max = 20L,
  inner_length_min = 18L, inner_length_max = 20L,
  tm_min = 55, tm_max = 72, gc_min = 35, gc_max = 65,
  dimer_dg_min = -8, hairpin_dg_min = -6, heterodimer_dg_min = -9,
  control_amplicon_min = 150L, control_amplicon_max = 200L,
  allele_amplicon_min = 80L, allele_amplicon_max = 150L,
  min_band_diff = 20L, max_candidates_per_pool = 8L,
  max_raw_candidates_per_pool = 8L
)
```

Ese caso tiene bandas esperadas de 152 bp (control), 105 bp (referencia) y
84 bp (alternativa). El límite de dímero cruzado es deliberadamente relajado;
revise el set con una herramienta termodinámica independiente antes de usarlo.

## Productos, gel y mapas

```r
pcr <- simulate_arms_pcr(
  gene_seq, design$best_set,
  snp_pos = snp$snp_pos,
  ref_allele = snp$snp_base,
  alt_allele = "A"
)

pcr$products
simulate_arms_gel(pcr)
export_arms_primers_txt(design, "arms_design.txt")
export_arms_amplicon_map_txt(pcr, "arms_nucleotide_map.txt")
```

`ref/ref` muestra control + referencia; `ref/alt`, las tres bandas; y
`alt/alt`, control + alternativa. El mapa TXT lista la hebra forward 5'→3' y
la complementaria 3'→5' alineadas base por base, con las coordenadas del
FASTA y los primers almacenados 5'→3'. Para el control heterocigoto se indica
que la banda representa ambas plantillas; la secuencia impresa usa la
representación de referencia.

## Pipeline

Cuando el flujo paso a paso esté revisado, use:

```r
result <- run_arms_pcr_pipeline(
  gene_seq,
  flank_seq = "SECUENCIA_5_PRIMA_ANTES_DEL_SNP",
  alt_allele = "A",
  output_file = "arms_design.txt"
)
```

El objeto conserva `snp`, `design`, `report`, `pcr_result` y `gel` para una
revisión reproducible.
