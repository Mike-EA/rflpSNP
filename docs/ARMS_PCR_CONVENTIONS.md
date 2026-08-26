# Convenciones de coordenadas y orientación para ARMS-PCR

Este documento fija el contrato de representación que usarán las funciones
ARMS-PCR. Es complementario al plan de desarrollo y no implementa el diseño de
cebadores.

## Secuencia y coordenadas

- `gene_seq` es un `Biostrings::DNAString` de la secuencia de referencia,
  leída y representada en dirección 5'→3'.
- Toda coordenada genómica (`snp_pos`, inicio y fin de un cebador) es
  **1-based**, inclusiva y se refiere a `gene_seq`, incluso si el flanco o el
  cebador se encuentra en la cadena complementaria.
- Un intervalo de `start` a `end` tiene longitud `end - start + 1`, con
  `start <= end`.
- `locate_snp()` devuelve `snp_pos` en este sistema. El campo `strand` explica
  dónde coincidió el flanco, pero no cambia el marco de coordenadas devuelto.

## Alelos

- `ref_allele` será la base en mayúscula de `gene_seq[snp_pos]` y deberá
  validarse contra esa base antes de diseñar un ensayo.
- `alt_allele` será un argumento explícito, una base `A`, `C`, `G` o `T`, en
  mayúscula y distinta de `ref_allele`. No se inferirá de la secuencia de
  referencia ni del flanco.

## Cebadores

- Todas las secuencias de cebadores se almacenarán y exportarán 5'→3'.
- Sus coordenadas describirán el segmento complementario en `gene_seq`:
  la orientación no invierte `start` y `end`.
- Un cebador de orientación `forward` coincide con
  `gene_seq[start:end]`. Un cebador `reverse` es
  `reverseComplement(gene_seq[start:end])`.
- Para un interno específico de alelo, el SNP coincide con el extremo 3':
  - el interno `forward` termina en `ref_allele` o `alt_allele` y su
    `end == snp_pos`;
  - el interno `reverse` termina en el complemento de `ref_allele` o
    `alt_allele` y su `start == snp_pos`.
- La discrepancia deliberada se identificará respecto al extremo 3' del
  cebador: `3'-2` y `3'-3` corresponden a los índices R `length(primer)-2` y
  `length(primer)-3`, respectivamente. El nucleótido terminal del SNP se
  denomina `3'` y no se cuenta como una discrepancia deliberada.

Estas reglas se comprobarán antes de construir candidatos en la etapa 1 y
antes de simular productos en la etapa 3.
