## Tips and tricks, lessons learned,...


1. n1 is not the right way to limit the number of multipoles

1. Scheme=3 is usually the most efficient and just as robust as Scheme=0

1. Balancing is quite useful, and has been extended to other schemes beyond Stout's

1. Masking blocks vs rows vs columns: unclear what is most meaningful at this stage, so all 3 are possible.

1. Cross-sections from collapsed collective. Seems to be more physically-interpretable in terms of individual orders, though it still produces negative partial cross-sections (except for scattering). Is this correct? The current interpretation leans towards a kind of interference between orders, in which case this is probably OK.

1. external T-matrices require the same wavelengths as in mode=2

1. TF, DF, etc are currently limited to 9

