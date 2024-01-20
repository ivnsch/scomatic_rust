# scomatic_rust

Port experiment of https://github.com/cortes-ciriano-lab/SComatic

[Step 1:](https://github.com/cortes-ciriano-lab/SComatic/blob/0f799b949ff1f32641bd6b7a3dde9a0ceb886aec/docs/SComaticExample.md#step-1-splitting-alignment-file-in-cell-type-specific-bams)

```
cargo run --release -- --bam ./test/Example.scrnaseq.bam \
    --meta ./test/Example.cell_barcode_annotations.tsv \
    --id Example \
    --n_trim 5 \
    --max_nM 5 \
    --max_NH 1 \
    --outdir .
```
