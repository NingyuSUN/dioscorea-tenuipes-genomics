# Pipeline figures

The figures describe analysis purposes and data relationships. They contain no sample-level results or claims that every module is automated or scientifically validated.

- `pipeline_overview.svg`: assembly, population and organelle branches, with later research modules.
- `pipeline_population.svg`: input audit, genotype panels, configured analyses and later sensitivity modules. Genotype-likelihood work uses source VCF PL fields and a matched SNP panel; PL likelihoods are not posterior probabilities.
- `pipeline_comparative.svg`: separate nuclear and chloroplast evidence streams, matched comparison and sensitivity checks, followed by conditional ILS scenarios. These scenarios are not a calibrated probability of a causal mechanism.

Line styles are part of the figure meaning:

- Solid teal: existing configured workflow connections.
- Dashed amber: separate historical scripts or templates.
- Dotted blue: research data relationships only; execution status and input contracts remain in the workflow documentation.

The overview is a method map, not a complete executable dependency graph. In particular, the later-research cards are separate modules, not a sequential chain. The comparative figure uses only dotted research connections.

Regenerate the three SVG files with Python 3 and its standard library:

```bash
python docs/figures/build_pipeline_figures.py
```

Edit labels, geometry and colors in `build_pipeline_figures.py`. Text remains editable SVG text, with 16 px as the smallest size. The files use no external images, scripts, fonts or web resources. Check the current workflow rules and `workflow/research/SOURCE_MAP.json` when changing module status or data inputs.
