# Letters

Short research notes, set in the preprint stack of the surrounding research
programme.

Each note here is the typeset form of an article published at
[gancor.xyz](https://gancor.xyz), which remains the canonical text. This
repository is the citable, archival copy: a PDF that can be cited and read
without depending on the site staying up, and, where a note has one, the code
that produced its numbers.

## Contents

| Note | Status | Date | Pages |
|---|---|---|---|
| [Loss-Versus-Rebalancing Under Emissions](letters/2026_05_17_lvr-under-emissions/lvr-under-emissions.pdf) | Technical Note | 17 May 2026 | 7 |

## Layout

One folder per note, named for its publication date so the series reads in
order.

```
letters/
  YYYY_MM_DD_<slug>/
    <slug>.pdf           the note
    figures/             figures it carries
    sim/                 verification code, where the note has any
    sim/OUTPUT.md        captured run, a regression target
```

A note without a `sim/` has nothing to reproduce.

## Licence

The MIT terms in `LICENSE` cover the code, the `sim/` folders.

The note PDFs and the figures under each note's `figures/` are (c) 2026
K. R. Ryan, all rights reserved. They may be read and cited; they may not be
redistributed or adapted without permission.

## Citing

Cite a note by its title and date, and link the article at gancor.xyz as the
canonical text. `CITATION.cff` describes the repository itself, for anyone
citing the code or the archival copy.
