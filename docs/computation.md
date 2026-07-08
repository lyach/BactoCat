# Computation of $k_{app}$ across enzyme cases

This document records how the apparent turnover number ($k_{app}$) are computed for each enzyme class, equations, units, and assumptions made. It follows Davidi et al. (2016).

A genome-scale metabolic model uses Gene-Protein-Reaction (GPR) rules to define which genes are required to produce the enzymes that catalyze a specific reaction. These are defined by `OR` and `AND` operators:

- ``AND-group``: represents a protein complex (heteromeric) where all specified subunits must be present together to function (e.g., b1 AND b2).

- ``OR-group``: represents isozymes, meaning they are alternative, standalone enzymes or complexes that can catalyze the exact same reaction independently (e.g., Enzyme A OR Enzyme B).

For each reaction in the model, its GPR rule is parsed into its **disjunctive
normal form (DNF)** using `gpr_to_dnf` in [`src/enzyme_classifier.py`](../src/enzyme_classifier.py). Each reaction is assigned an `enzyme_class`:

| `enzyme_class` | GPR shape | DNF example | Meaning |
|---|---|---|---|
| `homomeric` | single gene | `b1` -> `[{b1}]` | one protein catalyses the reaction |
| `complex` | one AND-group | `b1 and b2` -> `[{b1,b2}]` | one heteromeric holoenzyme |
| `isoenzyme` | OR of single genes | `b1 or b2` -> `[{b1},{b2}]` | interchangeable alternative enzymes |
| `mixed` | OR of AND-groups | `(b1 and b2) or b3` -> `[{b1,b2},{b3}]` | isozymes where at least one is a complex |

## Common quantities and units

For a gene $g$:

- $MW_g$ — molecular weight from the protein sequence, in g/mol (or mg/mmol). 
- $a_g$ — molar abundance `protein_mmol_gdcw`, in mmol/gDCW.
- $v$ — metabolic flux for the reaction, in mmol/gDCW/h (absolute value; reversibility/direction is handled upstream).

Unit conversions used:

- $v/3600$ converts mmol/gDCW/h to mmol/gDCW/s.
- $v\times 1000/60$ converts mmol/gDCW/h to umol/gDCW/min (for specific activity).

## Homomeric

A single gene catalyses the reaction. The apparent turnover number is the
 ratio of flux to molar enzyme abundance:

$$ k_{app}\,[\mathrm{s^{-1}}] = \frac{v/3600}{a_g} $$



## Isoenzymes and complexes

For multi-subunit (complex) reactions, the number of active sites is generally unknown. For multi-gene (isoenzyme) reactions, it is unlnkown how flux is partitioned. Following Davidi et al. (2016), we can treat the **whole
GPR as one integrated enzyme** and use the **enzyme mass fraction** as the
denominator.

### Mass fraction (denominator)

$$ M\,[\mathrm{mg/gDCW}] = \sum_{g} MW_g \cdot a_g $$

summed over the contributing polypeptide chains $g$. Because $MW_g$ is in mg/mmol and $a_g$ in
mmol/gDCW, the product is mg/gDCW.

- **Complex** (`b1 and b2 ...`): the sum runs over all subunits.
- **Isoenzyme** (`b1 or b2 ...`): the sum runs over all isozymes. The
  result is an "average enzyme" with a weighted kinetic constant.
- **Mixed** (`(b1 and b2) or b3`): the sum runs over every chain in every kept
  OR-branch, i.e. an isoenzyme that is also a complex contributes the sum of its subunit chains.

In all three cases the rule is simply: **sum $MW_g \cdot a_g$ over all genes in
the GPR** (subject to the missing-data policy outlined below).

### Specific activity

$$ SA_{app}\,[\mathrm{\mu mol\,mg^{-1}\,min^{-1}}] = \frac{v\times 1000/60}{M} $$

SA needs **no** assumption about active sites or stoichiometry.

### Turnover number ($s^{-1}$) conversion

To convert SA into a turnover number in $s^{-1}$, the value is computed directly from flux and the mass fraction:

$$ k_{app}\,[\mathrm{s^{-1}}] = \frac{(v/3600)\cdot MW_{tot}}{M}, \qquad MW_{tot} = \sum_g MW_g $$

where $MW_{tot}$ is the total enzyme weight (sum of contributing chain
weights).

Both `SA_app` and the converted `kcat_app` are reported (see the output columns
in [`docs/pipeline_specs.md`](pipeline_specs.md)).

## Assumptions

1. **Active sites = 1.** The conversion from specific activity to a per-enzyme
   turnover number divides the molecular weight by the number of catalytic
   sites. This value is unavailable for most complexes, so we assume one active
   site per complex. SA itself is unaffected; only the $s^{-1}$ value
   depends on this. 
2. **Subunit stoichiometry = 1.** Each polypeptide chain is counted exactly once in the sums for $M$ and $MW_{tot}$, without multiplication by copy number, because subunit stoichiometries are also unavailable.
3. **Integrated enzyme per reaction.** Because it is unknown how to partition flux among coexpressed isoenzymes (they may have different kinetics), they are merged into a single average enzyme. 
4. **Reporting both SA and $s^{-1}$.** SA is the assumption-free quantity; the $s^{-1}$ value is provided so isoenzymes/complexes sit on the same scale as homomeric enzymes for downstream comparison and ec-model use.

## Missing-data policy

Not every gene has a measured abundance ($a_g$) or a resolvable molecular weight
($MW_g$). A gene "has data" when both are present. The policy
(`_select_contributing_genes` in [`src/kapp_builder.py`](../src/kapp_builder.py))
differs by class:

- **Complex** (single AND-group): require **all** subunits to have data; if any
  subunit is missing, the reaction is dropped. A complex needs every
  subunit, so an incomplete sum would overestimate $SA$/$k_{app}$.
- **Isoenzyme** (OR of single genes): include **every** alternative that has
  data; an unmeasured isoenzyme contributes nothing. Drop
  only if no alternative has data.
- **Mixed** (OR of AND-groups): keep only the OR-branches whose subunits **all** have data (each such branch is a complete complex), then sum the chains of the kept branches. Drop if no branch is complete.

Enzymes with a non-positive or undefined mass fraction after applying the policy are dropped.

## $k_{max}$ and $\eta$

Identity for grouping across conditions is `enzyme_key`: the protein sequence for
homomeric enzymes and the integrated-enzyme id (`enzyme_id` = reaction + sorted
gene set) for the grouped classes.

- $k_{max} = \max_C k_{app}(C)$ — the maximum turnover across conditions
  (`kcat_app_max`), with the source condition recorded.
- $SA_{max} = \max_C SA_{app}(C)$ — reported as `SA_app_max`. For homomeric and
  fixed-composition enzymes this is the same condition as $k_{max}$ (since
  $MW_{tot}$ is constant), but for isozymes whose contributing set varies between
  conditions the two maxima may come from different conditions.
- $\eta(C) = k_{app}(C)/k_{max}$, with per-enzyme mean, stdev, min, max and
  coefficient of variation across conditions.

## Filtering

Apparent rates outside the physical range `[lower_threshold, upper_threshold]`
(default $10^{-5}$ to $10^{6}\ \mathrm{s^{-1}}$, the diffusion limit) are
discarded, for all enzyme classes, on the s$^{-1}$ turnover number.
